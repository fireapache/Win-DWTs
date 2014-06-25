#include <iostream>

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

void haar_Compress(double *vec, int n, double percentage)
{
	double energy = 0.0;
	double threshold;

	for (int i = 0; i < n; i++)
		energy += abs(vec[i]);

	threshold = energy * percentage;

	for (int i = 0; i < n; i++)
	if (abs(vec[i]) < threshold)
		vec[i] = 0.0;
}

#define POINTS 1024
#define PI atan(1) * 4

int main(int argc, char **argv)
{
	double vec[POINTS];
	double x1 = 0.0, x2 = 1.0;
	double x, alpha;

	for (int i = 0; i < POINTS; i++)
	{
		alpha = i / (double)(POINTS);
		x = x1 * (1.0 - alpha) + x2 * alpha;
		vec[i] = sin(2 * PI * x) + 1 - sqrt(abs(8 * PI * (x - 0.5)));
	}

	Haar_Decomposition(vec, 4, true);
	haar_Compress(vec, 4, 0.10);
	Haar_Composition(vec, 4, true);

	return 0;
}
