#pragma once
//-------------------------------------------------------------------------------------------------------------------------------
// This code comes from https://codes-sources.commentcamarche.net/source/27793-integration-numerique-methode-de-gauss-legendre
//-------------------------------------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <functional>

//double f(double x) /*fonction à intégrer */
/*{
	return(exp(-x * x / 2.0));
}*/

double pol_leg(int k, double x)
{
	int i;
	double p, p1, p2;
	p = 1; p1 = 0; p2 = 0;
	for (i = 1; i <= k; i++)
	{
		p2 = p1; p1 = p;
		p = ((2 * i - 1)*x*p1 - (i - 1)*p2) / i;
	}
	return(p);
}

double gauss_legendre(double a, double b, int k, std::function<double(double)> f)
{
	int j;
	double h, r, alpha, pol, s, t, racine[5];
	switch (k) {
	case 0:racine[0] = 0.0; break;
	case 1:racine[0] = 1.0 / sqrt(3.0); racine[1] = -racine[0]; break;
	case 2:racine[0] = sqrt(3.0 / 5.0); racine[1] = 0.0; racine[2] = -racine[0]; break;
	case 3:racine[0] = sqrt((15 + 2 * sqrt(30)) / 35); racine[1] = sqrt((15 - 2 * sqrt(30)) / 35);
		racine[2] = -racine[1]; racine[3] = -racine[0]; break;
	case 4:racine[0] = sqrt((35 + 2 * sqrt(70)) / 63); racine[1] = sqrt((35 - 2 * sqrt(70)) / 63);
		racine[2] = 0.0; racine[3] = -racine[1]; racine[4] = -racine[0]; break;
	}
	h = b - a;
	t = 0;
	s = 0.0;
	for (j = 0; j <= k; j++)
	{
		r = racine[j];
		pol = pol_leg(k, r);
		alpha = 2.0*(1.0 - r * r) / ((k + 1)*(k + 1)*pol*pol);
		s += alpha * f(h / 2 * (1 + r) + a);
	}
	t += h / 2 * s;
	return(t);
}
/*int main(int argc, char *argv[])
{
	int k;
	double a, b, t;

	a = 1; b = 2.0; // Borne d'intégration [a;b]

	printf("Méthode de Gauss-Legendre\n");

	printf("\nRang k      Valeur calculee       \n");
	for (k = 0; k <= 4; k++)
	{
		t = gauss_legendre(a, b, k);
		printf("%4d            %lf  \n", k, t);
	}

	return 0;
}*/

/*double integral(double a, double b, double(*f)(double))
{
	return gauss_legendre(a, b, 4, f);
}*/

double integral(double a, double b, std::function<double(double)> f)
{
	return gauss_legendre(a, b, 4, f);
}