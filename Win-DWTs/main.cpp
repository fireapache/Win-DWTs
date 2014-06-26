#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

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

#define ERROR 0.0000000001

void Haar_Compress(double *vec, int n, double percentage)
{
	double t, tMax, tMin, s, e;

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

	} while ((tMax - tMin) > ERROR);

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
	int levels = (int)log2((double)n);
	double *v;

	for (int i = n / 2; i > 1; i /= 2)
	{
		v = &vec[i - 1];
		Haar_Compress(v, i, percentage);
	}
}

void gnuplot_dat(const char *filename, double *x, double *y, int n)
{
	ofstream out;

	out.open(filename, ios_base::trunc);

	for (int i = 0; i < n; i++)
		out << x[i] << '\t' << y[i] << '\n';

	out.close();
}

#define POINTS 32
#define PI atan(1) * 4

int main(int argc, char **argv)
{
	double domain[POINTS];
	double image[POINTS];
	double x1 = 0.0, x2 = 1.0;
	double x, alpha;

	for (int i = 0; i < POINTS; i++)
	{
		alpha = i / (double)(POINTS);
		x = x1 * (1.0 - alpha) + x2 * alpha;
		domain[i] = x;
		image[i] = sin(2 * PI * x) + 1 - sqrt(abs(8 * PI * (x - 0.5)));
	}

	gnuplot_dat("test1.dat", domain, image, POINTS);

	Haar_Decomposition(image, POINTS, false);
	Haar_Levels_Compress(image, POINTS, 1.0);
	//Haar_Compress(image, POINTS, 0.90);
	Haar_Composition(image, POINTS, false);

	gnuplot_dat("test2.dat", domain, image, POINTS);

	return 0;
}
