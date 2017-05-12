// double evalNewton(double x, double* coeff, int n);
// double evalDeriv(double x, double* coeff, int n);
// double evalLegendre(double x, double* coeff, int n);
#include "new_legendreZeros.h"

using std::cout;
using std::setprecision;
using std::endl;

int getLegendre(double* zero, double* weights, double* coeff, int n)
{
	int i, count;
	if (errorCheck(coeff, n)) return 0;

	getLegendreCoeff(coeff,n);
	double TOL = 1E-16;

	for (i = 1; i < n+1; i++)
	{
		zero[i-1] = (1 - (1/(8*n*n)) + 1/(8*n*n*n)) * cos(M_PI*(4*i - 1)/(4*n+2));
	}

	for (i = 0; i < n; i++)
	{
		count = 0;
		double approx = zero[i];
		double error = 0;
		double func = 1;
		double deriv = 0;

		evaluate(approx,func,deriv,n);

		do
		{
			error = func / deriv;
			approx = approx - error;
			evaluate(approx,func,deriv,n);
			count ++;

		} while (fabs(error) > TOL && count < 10000);

		zero[i] = approx;
		weights[i] = 2/((1 - approx * approx) * deriv * deriv);
	}

	return 1;
}

void evaluate(double &approx, double &func, double &deriv, int n)
{
	double sub1 = approx;
	double sub2 = 1;
	double f = 1/ (approx * approx - 1);

	for(int i = 2; i <= n; ++i)
	{
		func = ((2 * i - 1) * approx * sub1 - (i - 1) * sub2) / ((double)i);
		deriv = i * f * (approx * func - sub1);

		sub2 = sub1;
		sub1 = func;
	}
}

/*
double evalNewton(double x, double* coeff, int n)
{
	double f = evalLegendre(x, coeff, n);
	double f_prime = evalDeriv(x, coeff, n);

	//printf("%.15f\n", f);
	//printf("%.15f\n", f_prime);
	//printf("%.15f\n",(f/f_prime));
	if (f_prime == 0.0) return x;
	return x + (double)(-1.0) * f/f_prime;
}

double evalLegendre(double x, double* coeff, int n)
{
	int i;
	double sum = 0;
	double xpwr = 1;

	for (i = 0; i < n+1; i++)
	{
		//printf("%.15f\n", sum);
		//printf("%.15f\n", coeff[i]);
		sum += xpwr * coeff[i];
		xpwr *= x;
	}

	//printf("%.15f\n", sum);
	return sum;
}

double evalDeriv(double x, double* coeff, int n)
{
	int i;
	double sum = 0;
	double xpwr = 1;

	for (i = 1; i < n+1; i++)
	{
		//printf("%.15f\n", coeff[i]);
		sum += xpwr * coeff[i] * (double)(i);
		xpwr *= x;
	}

	//printf("%.15f\n", sum);
	return sum;
}

*/

int getLegendreCoeff(double* coeff, int n)
{
	//coeff[n][i]
	// n = power
	// i = element
	n = n+1;

	double c[n+1][n+1];

	int cur;
	int i;

	for(cur = 0; cur < n+1; cur++)
	{
		for(i = 0; i < n+1; i++)
		{
			c[cur][i] = 0;
		}
	}

	c[0][0] = 1;

	if (n == 2)
	{
		coeff[0] = 1;
		return 1;
	}

	c[1][0] = 0;
	c[1][1] = 1;
	
	if (n == 1)
	{
		coeff[0] = 0;
		coeff[1] = 1;
		return 1;
	}

	for (cur = 2; cur < n; cur++)
	{
		for(i = 0; i < n; i++)
		{
			if (i-1 < 0)
				c[cur][i] = -c[cur-i][i]/cur;

			if (i > cur)
				c[cur][i] = 0;
			else
				c[cur][i] = ((2*cur - 1)*c[cur-1][i-1] - (cur-1)*c[cur-2][i])/cur;
		}
	}

	for (i = 0; i < n; i++)
	{
		coeff[i] = c[n-1][i];
	}

	return 1;
}

inline int errorCheck(double* coeff, int n)
{
	if (n <= 0)
	{
		printf("Parameter n out of bounds\n");
		return 1;
	}

	if (coeff == NULL)
	{
		printf("You must construct additional arrays!\n");
		return 1;
	}

	return 0;
}