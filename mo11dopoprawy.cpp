//:)

#include "stdafx.h"

#include<iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;
#pragma warning(disable:4996)

#define TOL 0.0000001            //dla metody iteracyjnej LU 10 ^ -6
#define TOLF 0.0000001
int loop = 50;
double tmax = 0.5;
double dt = 0.001;
double h = 0.1;
double lambda = dt / (h*h);

double rozw_analityczne(double x, double t)
{
	return 1.0 + exp(-1.0 * pow(M_PI, 2.0)*t)*cos(M_PI*x);
}

double war_poczatkowy(double x)
{
	return 1.0 + cos(M_PI*x);
}

double est(double *X, double *x_nowe)
{
	double x[4];
	x[0] = fabs(X[0] - x_nowe[0]);
	x[1] = fabs(X[1] - x_nowe[1]);
	x[2] = fabs(X[2] - x_nowe[2]);
	x[3] = fabs(X[3] - x_nowe[3]);

	double max = x[0];
	if (x[1] > max) max = x[1];
	if (x[2] > max) max = x[2];
	if (x[3] > max) max = x[3];

	return max;
}

double res(double * l, double *d, double *u, double *b, double *x_nowe, int rozmiar) //przekazuje macierz A i x nowe wyliczam wektor Ax-b i licze z niego maximum czyli najwieksza bezwzgledna wartosc  tego wektora
{
	double *Ax = new double[rozmiar];

	Ax[0] = fabs((d[0] * x_nowe[0] + u[0] * x_nowe[1]) - b[0]);
	for (int i = 1; i <= rozmiar - 2; i++)
	{
		Ax[i] = fabs((l[i - 1] * x_nowe[i - 1] + d[i] * x_nowe[i] + u[i] * x_nowe[i + 1]) - b[i]);
	}
	Ax[rozmiar - 1] = fabs((l[rozmiar - 2] * x_nowe[rozmiar - 2] + d[rozmiar - 1] * x_nowe[rozmiar - 1]) - b[rozmiar - 1]);

	double max = Ax[0];
	for (int i = 1; i < rozmiar; i++)
	{
		if (Ax[i] > max) max = Ax[i];
	}

	return max;
}

double resM(double * M, double *b, double *x_nowe, int rozmiar)
{
	return 0.0;
}

double maxblad(double *BLAD_thomas, int rozmiar)
{
	double max = BLAD_thomas[1];  // omijam z warunkow brzegowych

	for (int i = 1; i < rozmiar; i++)   //szukam najwiekszego bledu
	{
		if (BLAD_thomas[i] > max) max = BLAD_thomas[i];
	}

	return max;
}

void uzupelnijLDUB(double * l, double *d, double *u, double *vecb, double *Uk, int rozmiar)
{
	for (int i = 0; i <rozmiar - 2; i++)
	{
		l[i] = lambda / 2.0;
	}
	l[rozmiar - 2] = -1.0 / h;

	u[0] = 1.0 / h;
	for (int i = 1; i < rozmiar - 1; i++)
	{
		u[i] = lambda / 2.0;
	}

	d[0] = -1.0 / h;
	for (int i = 1; i < rozmiar - 1; i++)
	{
		d[i] = -1.0*(1.0 + lambda);
	}
	d[rozmiar - 1] = 1.0 / h;

	int index = 0;
	vecb[0] = 0.0;

	for (int i = 1; i < rozmiar - 1; i++)
	{
		vecb[i] = (-lambda / 2.0)*Uk[index] - (1.0 - lambda)*Uk[index + 1] - 1.0 * (lambda / 2.0)*Uk[index + 2];
		index++;
	}
	vecb[rozmiar - 1] = 0.0;

}

void uzupb(double * l, double *d, double *u, double *vecb, double *Uk, int rozmiar)
{
	vecb[0] = 0.0;
	int index = 0;
	for (int i = 1; i < rozmiar - 1; i++)
	{
		vecb[i] = (-lambda / 2.)*Uk[index] - (1.0 - lambda)*Uk[index + 1] - 1.0 * (lambda / 2.0)*Uk[index + 2];
		index++;
	}
	vecb[rozmiar - 1] = 0.0;
}

double * AlgorytmThomasa(double * l, double *d, double *u, double *b, double *x, int rozmiar)
{
	double *ni = new double[rozmiar];
	double *bb = new double[rozmiar];
	ni[0] = d[0];

	for (int i = 0; i < rozmiar - 1; i++)
	{
		ni[i + 1] = d[i + 1] - ((l[i] * u[i]) / ni[i]);

	}

	bb[0] = b[0];
	for (int i = 0; i < rozmiar - 1; i++)
		bb[i + 1] = b[i + 1] - ((l[i] * bb[i]) / ni[i]);

	x[rozmiar - 1] = bb[rozmiar - 1] / ni[rozmiar - 1];

	for (int i = rozmiar - 2; i >= 0; i--)
	{
		x[i] = (bb[i] - u[i] * x[i + 1]) / ni[i];
	}

	delete[] ni;
	delete[] bb;

	return x;
}

void swapRows(double **matrix, int j1, int j2, double *vector, int rozmiar)
{
	double tmp;

	for (int i = 0; i < rozmiar; i++)
	{
		tmp = matrix[j1][i];
		matrix[j1][i] = matrix[j2][i];
		matrix[j2][i] = tmp;
	}

	tmp = vector[j1];
	vector[j1] = vector[j2];
	vector[j2] = tmp;
}

int findAbsMax(double **matrix, int i, int j, int rozmiar)
{
	int k = 0;
	int indexMax = i;
	if (indexMax == rozmiar - 1)
	{
		return indexMax;
	}

	for (k = i + 1; k < rozmiar; k++)
	{
		if (abs(matrix[k][j]) > abs(matrix[indexMax][j]))
		{
			indexMax = k;
		}
	}

	return indexMax;
}

void partialChoice(double **matrix, int i, int j, double *vector, int rozmiar)
{
	swapRows(matrix, findAbsMax(matrix, i, j, rozmiar), j, vector, rozmiar);
}

void solveLowerTriangular(double **matrix, double *y, double *b, int rozmiar)
{
	for (int i = 0; i < rozmiar; i++)
	{
		double sum = 0.0;
		for (int j = 0; j <= i - 1; j++)
		{
			sum += matrix[i][j] * y[j];
		}
		y[i] = (b[i] - sum); // /1.0 (bo na przekatnej sa jedynki)
	}
}

void solveUpperTriangular(double **matrix, double *x, double *b, int rozmiar)
{
	for (int i = rozmiar - 1; i >= 0; i--)
	{
		double sum = 0.0;
		for (int j = i + 1; j < rozmiar; j++)
		{
			sum += matrix[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / matrix[i][i];
	}
}

void gaussianEliminate(double **matrix, int ip, int jp, int rozmiar)
{
	for (int i = ip + 1; i < rozmiar; i++)
	{
		matrix[i][jp] = matrix[i][jp] / matrix[ip][jp];
		for (int j = jp + 1; j < rozmiar; j++)
		{
			matrix[i][j] = matrix[i][j] - matrix[i][jp] * matrix[ip][j];
		}
	}

}

void printMatrix(double **matrix, int rozmiar)
{
	int i = 0, j = 0;

	for (i = 0; i < rozmiar; i++)
	{
		for (j = 0; j < rozmiar; j++)
		{
			printf("%f ", matrix[i][j]);
		}
		cout << "\n";
	}
	cout << endl << endl;
}

double * LU(double * l, double *d, double *u, double *b, int rozmiar, double **M)
{
	double sum = 0.0, *X;
	int counter = 0;
	int index = 0;
	double *xx = new double[rozmiar];
	double *yy = new double[rozmiar];

	for (int i = 0; i < rozmiar; i++)
	{
		if (M[i][i] == 0.0)
		{
			partialChoice(M, i, i, b, rozmiar);
		}

		gaussianEliminate(M, i, i, rozmiar);
	}

	X = new double[rozmiar];

	for (int i = 0; i < rozmiar; i++)//wypelniam wektor pryblizen poczatkowych
		X[i] = 0.0;

	while (1)
	{
		double estymator;
		double residuum;

		solveLowerTriangular(M, yy, b, rozmiar);//uzupelniam yy
		solveUpperTriangular(M, xx, yy, rozmiar);//uzupe³niam xx

		estymator = est(X, xx);    //licze estymator jako najwieksza roznica pomiedzy odpowiednimi elementami tych wektorow
		residuum = res(l, d, u, b, xx, rozmiar);

		for (int i = 0; i < rozmiar; i++)                    //przepisuje nowowyliczone X[i] z tymczasowego wektora newX[i]
		{
			X[i] = xx[i];
			if (fabs(X[i]) > 1000)
			{
				cout << 'a';
			}
		}

		if (counter >= loop) 
			break;

		if (estymator < TOL) 
			break;
		
		if (residuum < TOLF) 
			break;

		counter++;
	}

	delete[] X;
	return xx;
}

void rozwiaz_rownanie(double a, double b, double h)
{
	ofstream file("dane.csv", std::ofstream::trunc);
	locale loc = locale("");
	file.imbue(loc);//zmieniam locale pliku, zeby w excelu byl odpowiedni separator(import csv)
	if (!file.is_open())
	{
		cout << "Blad otwarcia pliku" << endl;
		return;
	}

	ofstream file2("blad.csv", std::ofstream::trunc);
	file2.imbue(loc);
	if (!file2.is_open())
	{
		cout << "Blad otwarcia pliku" << endl;
		return;
	}

	int rozmiar = static_cast<int>(fabs(a - b) / h) + 1;
	double *l = new double[rozmiar - 1];
	double *d = new double[rozmiar];
	double *u = new double[rozmiar - 1];
	double *vecb;
	double *Uk_thomas = new double[rozmiar];
	double *Uk_LU = new double[rozmiar];

	double *BLAD_thomas = new double[rozmiar];
	double *BLAD_LU = new double[rozmiar];

	int index = 0;

	//tworze macierze l,d,u i vektor b
	double **M = new double*[rozmiar];
	for (int j = 0; j<rozmiar; j++)
	{
		M[j] = new double[rozmiar];
	}

	vecb = new double[rozmiar + 1];
	double *x = new double[rozmiar + 1];

	printf("T\t NUMERYCzNIE \t ANALITYCZNIE");

	for (double x = a; x < b; x += h)    //uzupelniam Uk z war poczatkowego , potrzebne do vectora b na starcie programu
	{
		Uk_thomas[index] = war_poczatkowy(x);
		Uk_LU[index] = war_poczatkowy(x);

		index++;
	}

	uzupelnijLDUB(l, d, u, vecb, Uk_thomas, rozmiar);
	uzupelnijLDUB(l, d, u, vecb, Uk_LU, rozmiar);

	for (int k = 0; k < rozmiar; k++)
	{
		for (int i = 0; i < rozmiar; i++)
		{
			if (k == i + 1)
				M[k][i] = l[i];
			else
			{
				if (i == k + 1)
					M[k][i] = u[k];
				else
				{
					if (k == i)
						M[k][i] = d[k];
					else
						M[k][i] = 0.0;
				}
			}
		}
	}

	for (double t = 0.0; t <= tmax; t += dt)
	{
		uzupb(l, d, u, vecb, Uk_thomas, rozmiar);

		Uk_thomas = AlgorytmThomasa(l, d, u, vecb, x, rozmiar);
		Uk_LU = LU(l, d, u, vecb, rozmiar, M);

		printf("   T\t  h\t  X  \t    THOMAS \t   LU              ANALITYCZNIE:\tBLAD(thomas) \t BLAD(LU)\n");

		index = 0;
		for (double x = a; x < b; x += h)//wypisuje wyniki dla tego t
		{
			BLAD_thomas[index] = fabs(rozw_analityczne(x, t) - Uk_thomas[index]);
			BLAD_LU[index] = fabs(rozw_analityczne(x, t) - Uk_LU[index]);

			printf("%6.4f\t %6.4f\t %6.4f\t     %6.6f \t %6.6f \t %6.6f \t %6.6f \t %6.6f \n", t, h, x, Uk_thomas[index], Uk_LU[index], rozw_analityczne(x, t), BLAD_thomas[index], BLAD_LU[index]);

			//file << t << ";" << h << ";" << x << ";" << Uk_thomas[index] << ";" << Uk_LU[index] << ";" << rozw_analityczne(x, t) << ";" << BLAD_thomas[index] << ";" << BLAD_LU[index] << endl;
			index++;
		}

		file2 << t << ";" << maxblad(BLAD_thomas, rozmiar) << ";" << maxblad(BLAD_LU, rozmiar) << ";" << endl;

	}

	delete[] l, d, u, vecb, Uk_thomas;

	for (int i = 0; i < rozmiar; i++)
		delete[](M[i]);
	delete[](M);

	file.close();
	file2.close();
}


int main()
{
	rozwiaz_rownanie(0.0, 1.0, h);
	getchar();
	return 0;
}

