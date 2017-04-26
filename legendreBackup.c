#include <stdio.h>
#include <math.h>
#include <time.h>

int getLegendreCoeff(double* a, int n);
int getZeroes(double* zero, double* a, int n);
double evalNewton(double x, double* a, int n);
double evalDeriv(double x, double* a, int n);
double evalLegendre(double x, double* a, int n);
int errorCheck(double* a, int n);

int getZeroes(double* zero, double* a, int n)
{
	int i, count;
	if (errorCheck(a, n)) return 0;

	//printf("before LegCoeff\n");
	//for (i = 0; i<n; i++)printf("%.15f\n", a[i]);
	//printf("end LegCoeff\n");
	getLegendreCoeff(a,n);

	double TOL = 1E-16;
	//printf("after legCoeff\n");
	//for (i = 0; i < n+1; i++)	printf("%.15f\n", a[i]);
	//printf("end legCoeff\n");

	for (i = 1; i < n+1; i++)
	{
		zero[i-1] = (1 - (1/(8*n*n)) + 1/(8*n*n*n)) * cos(M_PI*(4*i - 1)/(4*n+2));
		//printf("%.15f\n", zero[i-1]);
	}

	for (i = 0; i < n; i++)
	{
		count = 0;
		double approx = zero[i];
		//printf("%.15f\n", zero[i]);
		double check;
		double error;

		do
		{
			//printf("%.15f\n", approx);
			check = evalNewton(approx, a, n);
			//printf("%.15f\n", check);
			error = fabs(check - approx);
			//printf("%.15f\n", error);
			approx = check;
			//printf("%.15f\n", check);
			count ++;

		} while (error > TOL && count < 1000000000);

		zero[i] = approx;
		//printf("%.15f\n", zero[i]);
	}

	//printf("end legZeros\n");
	for (i = 0; i < n; i++)
		printf("The %ith root is %.15f\n", i, zero[i]);

	return 1;
}


double evalNewton(double x, double* a, int n)
{
	double f = evalLegendre(x, a, n);
	double f_prime = evalDeriv(x, a, n);

	//printf("%.15f\n", f);
	//printf("%.15f\n", f_prime);
	//printf("%.15f\n",(f/f_prime));
	if (f_prime == 0.0) return x;
	return x + (double)(-1.0) * f/f_prime;
}

double evalLegendre(double x, double* a, int n)
{
	int i;
	double sum = 0;
	double xpwr = 1;

	for (i = 0; i < n+1; i++)
	{
		//printf("%.15f\n", sum);
		//printf("%.15f\n", a[i]);
		sum += xpwr * a[i];
		xpwr *= x;
	}

	//printf("%.15f\n", sum);
	return sum;
}

double evalDeriv(double x, double* a, int n)
{
	int i;
	double sum = 0;
	double xpwr = 1;

	for (i = 1; i < n+1; i++)
	{
		//printf("%.15f\n", a[i]);
		sum += xpwr * a[i] * (double)(i);
		xpwr *= x;
	}

	//printf("%.15f\n", sum);
	return sum;
}



int getLegendreCoeff(double* a, int n)
{
	//a[n][i]
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
		a[0] = 1;
		return 1;
	}

	c[1][0] = 0;
	c[1][1] = 1;
	
	if (n == 1)
	{
		a[0] = 0;
		a[1] = 1;
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
		a[i] = c[n-1][i];
	}

	return 1;
}

inline int errorCheck(double* a, int n)
{
	if (n <= 0)
	{
		printf("Parameter n out of bounds\n");
		return 1;
	}

	if (a == NULL)
	{
		printf("You must construct additional arrays!\n");
		return 1;
	}

	return 0;
}