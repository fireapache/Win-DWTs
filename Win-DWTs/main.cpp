#include <iostream>
#include <cmath>
#include <fstream>
#include <windows.h>

using namespace std;

double PCFreq = 0.0;
__int64 CounterStart = 0;

void StartCounter()
{
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li))
		cout << "QueryPerformanceFrequency failed!\n";

	PCFreq = double(li.QuadPart) / 1000.0;

	QueryPerformanceCounter(&li);
	CounterStart = li.QuadPart;
}

double GetCounter()
{
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return double(li.QuadPart - CounterStart) / PCFreq;
}

void ResetCounter()
{
	PCFreq = 0.0;
	CounterStart = 0;
}

// Passo de composição tipo Haar de um vetor.
void Haar_CompositionStep(double *vec, int n, bool normal)
{
	double *vecp = new double[2 * n];

	for (int i = 0; i < n; i++)
		vecp[i] = 0;

	for (int i = 0; i < n; i++)
	{
		vecp[2 * i] = vec[i] + vec[n + i];
		vecp[2 * i + 1] = vec[i] - vec[n + i];

		if (normal)
		{
			vecp[2 * i] /= sqrt(2.0);
			vecp[2 * i + 1] /= sqrt(2.0);
		}
	}

	for (int i = 0; i < 2 * n; i++)
		vec[i] = vecp[i];

	delete[] vecp;
}

void VinisNormalization(double *vec, UINT n)
{
	UINT levels = log2(n);

	for (UINT level = 0; level < levels; level++)
	{
		UINT inicio;

		if (level == 0) inicio = 0;
		else        	inicio = (UINT)pow(2.0, (double)(level));

		UINT fim = (UINT)pow(2.0, (double)(level + 1));

		for (UINT i = inicio; i < fim; i++)
			vec[i] *= pow(2.0, -(double)level / 2.0);
	}
}

void VinisStandardMatrixNormalization(double **mat, UINT n, bool invert)
{
	UINT levels = log2(n);

	for (UINT levelL = 0; levelL < levels; levelL++)
	{
		UINT inicioL;

		if (levelL == 0)	inicioL = 0;
		else        		inicioL = (UINT)pow(2.0, (double)(levelL));

		UINT fimL = (UINT)pow(2.0, (double)(levelL + 1));

		for (UINT l = inicioL; l < fimL; l++)
		{
			for (UINT levelC = 0; levelC < levels; levelC++)
			{
				UINT inicioC;

				if (levelC == 0)	inicioC = 0;
				else        		inicioC = (UINT)pow(2.0, (double)(levelC));

				UINT fimC = (UINT)pow(2.0, (double)(levelC + 1));

				if (invert)
				{
					for (UINT c = inicioC; c < fimC; c++)
						mat[l][c] *= pow(2.0, (double)(levelL + levelC) / 2.0);
				}
				else
				{
					for (UINT c = inicioC; c < fimC; c++)
						mat[l][c] /= pow(2.0, (double)(levelL + levelC) / 2.0);
				}
			}
		}
	}
}

void VinisNonStandardMatrixNormalization(double **matrix, UINT n, bool invert)
{
	double div;
	UINT start, limit;

	UINT level = (UINT)(log((double)n) / log(2.0)) - 1;

	if (level <= 0) return;

	for (UINT i = level; i > 0; i--)
	{
		div = pow(2.0, (double)i);
		start = (UINT)div;
		limit = (UINT)pow(2.0, (double)i + 1);

		for (UINT l = 0; l < limit / 2; l++)
		{
			for (UINT c = start; c < limit; c++)
			{
				matrix[l][c] /= div;
			}
		}

		if (invert)
		{
			for (UINT l = start; l < limit; l++)
			for (UINT c = 0; c < limit; c++)
				matrix[l][c] *= div;
		}
		else
		{
			for (UINT l = start; l < limit; l++)
			for (UINT c = 0; c < limit; c++)
				matrix[l][c] /= div;
		}

	}
}

void VinisMatrixNormalization(double **mat, UINT n, bool standard, bool invert)
{
	if (standard)	VinisStandardMatrixNormalization(mat, n, invert);
	else			VinisNonStandardMatrixNormalization(mat, n, invert);
}

// Composição completa tipo Haar de um vetor.
void Haar_Composition(double *vec, int n, bool normal)
{
	for (int i = 1; i < n; i = i * 2)
		Haar_CompositionStep(vec, i, normal);

	if (normal)
	for (int i = 0; i < n; i++)
		vec[i] = vec[i] * sqrt(float(n));
}

void Haar_DecompositionStep(double *vec, int n, bool normal)
{
	double div;
	double *vecp = new double[n];

	for (int i = 0; i < n; i++)
		vecp[i] = 0;

	if (normal) div = sqrt(2.0);
	else        div = 2.0;

	for (int i = 0; i < n / 2; i++)
	{
		vecp[i] = (vec[2 * i] + vec[2 * i + 1]) / div;
		vecp[i + n / 2] = (vec[2 * i] - vec[2 * i + 1]) / div;
	}

	for (int i = 0; i < n; i++)
		vec[i] = vecp[i];

	delete[] vecp;
}

// Decomposição completa de tipo Haar de um vetor.
void Haar_Decomposition(double *vec, int n, bool normal)
{
	if (normal)
	for (int i = 0; i < n; i++)
		vec[i] = vec[i] / sqrt(float(n));

	while (n > 1)
	{
		Haar_DecompositionStep(vec, n, normal);
		n /= 2;
	}
}

#define THRESHOLD_ERROR 0.0000000001

void Haar_Compress(double *vec, int n, double percentage)
{
	double t, tMax, tMin, s, e;

	if (percentage >= 1.0)
	{
		for (int i = 0; i < n; i++)
		{
			vec[i] = 0.0;
		}
	}

	tMax = tMin = abs(vec[0]);
	s = 0.0;

	for (int i = 0; i < n; i++)
	{
		if (abs(vec[i]) > tMax) tMax = abs(vec[i]);
		if (abs(vec[i]) < tMin) tMin = abs(vec[i]);
		s += abs(vec[i]);
	}

	e = s * percentage;

	do
	{
		t = (tMax + tMin) / 2.0;
		s = 0.0;

		for (int i = 0; i < n; i++)
		{
			if (abs(vec[i]) < t) s += pow(vec[i], 2.0);
		}

		if (s < pow(e, 2.0)) tMin = t;
		else tMax = t;

	} while ((tMax - tMin) > THRESHOLD_ERROR);

	for (int i = 0; i < n; i++)
	{
		if (abs(vec[i]) < t)
		{
			vec[i] = 0.0;
		}
	}
}

void Haar_Levels_Compress(double *vec, int n, double percentage)
{
	double *v;

	for (int i = 1; i < n; i *= 2)
	{
		v = &vec[i];
		Haar_Compress(v, i, percentage);
	}
}

void Haar_StandardDecomposition(double **matrix, int rows, int cols, bool normal)
{
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
			temp_row[j] = matrix[i][j];

		Haar_Decomposition(temp_row, cols, normal);

		for (int j = 0; j < cols; j++)
			matrix[i][j] = temp_row[j];
	}

	for (int i = 0; i < cols; i++)
	{
		for (int j = 0; j < rows; j++)
			temp_col[j] = matrix[j][i];

		Haar_Decomposition(temp_col, rows, normal);

		for (int j = 0; j < rows; j++)
			matrix[j][i] = temp_col[j];
	}

	delete[] temp_row;
	delete[] temp_col;
}

void Haar_NonStandardDecomposition(double **matrix, int rows, int cols, bool normal)
{
	int h = rows, w = cols;
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];


	if (normal)
	{
		for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			matrix[i][j] = matrix[i][j] / cols;
	}


	while (w > 1 || h > 1)
	{
		if (w > 1)
		for (int i = 0; i < h; i++)
		{
			for (int j = 0; j < w; j++)
				temp_row[j] = matrix[i][j];

			Haar_DecompositionStep(temp_row, w, normal);

			for (int j = 0; j < w; j++)
				matrix[i][j] = temp_row[j];
		}

		if (h > 1)
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
				temp_col[j] = matrix[j][i];

			Haar_DecompositionStep(temp_col, h, normal);

			for (int j = 0; j < h; j++)
				matrix[j][i] = temp_col[j];
		}

		if (w > 1) w /= 2;
		if (h > 1) h /= 2;
	}
}

void Haar_MatrixDecomposition(double **matrix, int rows, int cols, bool normal, bool standard)
{
	/*
	if (rows != cols)
	{
	cout << "Não é uma matriz quadrada!" << endl;
	return;
	}
	*/

	if (standard) Haar_StandardDecomposition(matrix, rows, cols, normal);
	else          Haar_NonStandardDecomposition(matrix, rows, cols, normal);
}

void Haar_StandardComposition(double **matrix, int rows, int cols, bool normal)
{
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];

	for (int i = 0; i < cols; i++)
	{
		for (int j = 0; j < rows; j++)
			temp_col[j] = matrix[j][i];

		Haar_Composition(temp_col, rows, normal);

		for (int j = 0; j < rows; j++)
			matrix[j][i] = temp_col[j];
	}

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
			temp_row[j] = matrix[i][j];

		Haar_Composition(temp_row, cols, normal);

		for (int j = 0; j < cols; j++)
			matrix[i][j] = temp_row[j];
	}

	delete[] temp_row;
	delete[] temp_col;
}

void Haar_NonStandardComposition(double **matrix, int rows, int cols, bool normal)
{
	int h = 1, w = 1;
	double *temp_row = new double[cols];
	double *temp_col = new double[rows];

	while (w < cols || h < rows)
	{
		if (h < rows)
		for (int i = 0; i < 2 * w; i++)
		{
			for (int j = 0; j < rows; j++)
				temp_col[j] = matrix[j][i];

			Haar_CompositionStep(temp_col, h, normal);

			for (int j = 0; j < rows; j++)
				matrix[j][i] = temp_col[j];
		}

		if (w < cols)
		for (int i = 0; i < 2 * h; i++)
		{
			for (int j = 0; j < cols; j++)
				temp_row[j] = matrix[i][j];

			Haar_CompositionStep(temp_row, w, normal);

			for (int j = 0; j < cols; j++)
				matrix[i][j] = temp_row[j];
		}

		if (w < cols) w = w * 2;
		if (h < rows) h = h * 2;
	}


	if (normal)
	{
		for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			matrix[i][j] = matrix[i][j] * cols;
	}

	delete[] temp_row;
	delete[] temp_col;
}

void Haar_MatrixComposition(double **matrix, int rows, int cols, bool normal, bool standard)
{
	/*
	if (rows != cols)
	{
	cout << "Não é uma matriz quadrada!" << endl;
	return;
	}
	*/

	if (standard) Haar_StandardComposition(matrix, rows, cols, normal);
	else          Haar_NonStandardComposition(matrix, rows, cols, normal);
}

void gnuplot_dat(const char *filename, double *x, double *y, int n)
{
	ofstream out;

	out.open(filename, ios_base::trunc);

	for (int i = 0; i < n; i++)
		out << x[i] << '\t' << y[i] << '\n';

	out.close();
}

void gnuplot_dat_Vdecomposition(const char *file, double x1, double x2, double *v, int n, int levels, bool normal)
{
	ofstream out;
	double x, alpha;

	if (normal)
	for (int i = 0; i < n; i++)
		v[i] = v[i] / sqrt(float(n));

	for (int i = 0; i < levels; i++)
	{
		Haar_DecompositionStep(v, n, normal);
		n /= 2;
	}

	out.open(file, ios_base::trunc);

	for (int i = 0; i < n; i++)
	{
		alpha = (double)i / (double)(n);
		x = x1 * (1.0 - alpha) + x2 * alpha;

		out << x << '\t' << v[i] << '\n';
	}

	out.close();

	for (int i = levels; i > 0; i--)
	{
		Haar_CompositionStep(v, n, normal);
		n *= 2;
	}

	if (normal)
	for (int i = 0; i < n; i++)
		v[i] = v[i] * sqrt(float(n));
}

void gnuplot_dat_Wdecomposition(const char *file, double x1, double x2, double *v, int n, int levels, bool normal)
{
	ofstream out;
	double x, alpha;

	if (normal)
	for (int i = 0; i < n; i++)
		v[i] = v[i] / sqrt(float(n));

	for (int i = 0; i < levels; i++)
	{
		Haar_DecompositionStep(v, n, normal);
		n /= 2;
	}

	out.open(file, ios_base::trunc);

	for (int i = 0; i < n; i++)
	{
		alpha = (double)i / (double)(n);
		x = x1 * (1.0 - alpha) + x2 * alpha;
		x += (x2 - x1) / (2 * n);
		out << x << '\t' << v[i + n] << '\n';
	}

	out.close();

	for (int i = levels; i > 0; i--)
	{
		Haar_CompositionStep(v, n, normal);
		n *= 2;
	}

	if (normal)
	for (int i = 0; i < n; i++)
		v[i] = v[i] * sqrt(float(n));
}

void gnuplot_dat_VWdecomposition(const char *file1, const char *file2, double x1, double x2, double *v, int n, int levels, bool normal)
{
	ofstream out1, out2;
	double x, alpha;

	if (normal)
	for (int i = 0; i < n; i++)
		v[i] = v[i] / sqrt(float(n));

	for (int i = 0; i < levels; i++)
	{
		Haar_DecompositionStep(v, n, normal);
		n /= 2;
	}

	out1.open(file1, ios_base::trunc);
	out2.open(file2, ios_base::trunc);

	for (int i = 0; i < n; i++)
	{
		alpha = (double)i / (double)(n);
		x = x1 * (1.0 - alpha) + x2 * alpha;

		out1 << x << '\t' << v[i] << '\n';
		out2 << x << '\t' << v[i + n] << '\n';
	}

	out1.close();
	out2.close();

	for (int i = levels; i > 0; i--)
	{
		Haar_CompositionStep(v, n, normal);
		n *= 2;
	}

	if (normal)
	for (int i = 0; i < n; i++)
		v[i] = v[i] * sqrt(float(n));
}

void test0(); // 27/06/2014
void test1(); // 22/07/2014

int main(int argc, char **argv)
{
	test0();

	return 0;
}

#define POINTS 1024

void test1()
{

}

#define PI atan(1) * 4

void test0()
{
	double domain[POINTS];
	double image[POINTS];
	double t1[POINTS], t2[POINTS], t3[POINTS];
	double x1 = 0.0, x2 = 1.0;
	double x, alpha;

	for (int i = 0; i < POINTS; i++)
	{
		alpha = (double)i / (double)(POINTS);
		x = x1 * (1.0 - alpha) + x2 * alpha;
		domain[i] = x;
		image[i] = sin(2 * PI * x) + 1 - sqrt(abs(8 * PI * (x - 0.5)));
		t1[i] = t2[i] = t3[i] = image[i];
	}

	gnuplot_dat_Vdecomposition("1_V.dat", 0.0, 1.0, image, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("1_1W.dat", 0.0, 1.0, image, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("1_2W.dat", 0.0, 1.0, image, POINTS, 2, true);
	gnuplot_dat_Wdecomposition("1_3W.dat", 0.0, 1.0, image, POINTS, 3, true);

	Haar_Decomposition(t1, POINTS, true);
	//Haar_Compress(t1, POINTS, 0.01);
	Haar_Levels_Compress(t1, POINTS, 0.01);
	Haar_Composition(t1, POINTS, true);

	gnuplot_dat_Vdecomposition("2_V.dat", 0.0, 1.0, t1, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("2_1W.dat", 0.0, 1.0, t1, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("2_2W.dat", 0.0, 1.0, t1, POINTS, 2, true);
	gnuplot_dat_Wdecomposition("2_3W.dat", 0.0, 1.0, t1, POINTS, 3, true);

	Haar_Decomposition(t2, POINTS, true);
	//Haar_Compress(t2, POINTS, 0.02);
	Haar_Levels_Compress(t2, POINTS, 0.02);
	Haar_Composition(t2, POINTS, true);

	gnuplot_dat_Vdecomposition("3_V.dat", 0.0, 1.0, t2, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("3_1W.dat", 0.0, 1.0, t2, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("3_2W.dat", 0.0, 1.0, t2, POINTS, 2, true);
	gnuplot_dat_Wdecomposition("3_3W.dat", 0.0, 1.0, t2, POINTS, 3, true);

	Haar_Decomposition(t3, POINTS, true);
	//Haar_Compress(t3, POINTS, 0.035);
	Haar_Levels_Compress(t3, POINTS, 0.035);
	Haar_Composition(t3, POINTS, true);

	gnuplot_dat_Vdecomposition("4_V.dat", 0.0, 1.0, t3, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("4_1W.dat", 0.0, 1.0, t3, POINTS, 1, true);
	gnuplot_dat_Wdecomposition("4_2W.dat", 0.0, 1.0, t3, POINTS, 2, true);
	gnuplot_dat_Wdecomposition("4_3W.dat", 0.0, 1.0, t3, POINTS, 3, true);
}