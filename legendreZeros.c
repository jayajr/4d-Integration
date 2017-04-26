#include <stdio.h>
#include <math.h>

int getLegendreCoeff(double* a, int n);
int getZeroes(double* zero, double* a, int n);
int errorCheck(double* a, int n);

int getZeroes(double* zero, double* a, int n)
{
	int i,j, count;
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
		double deriv;
		double error;





		do
		{
			double acc1  = 0;		double acc11 = 0;
			double acc2  = 0;		double acc12 = 0;
			double acc3  = 0;		double acc13 = 0;
			double acc4  = 0;		double acc14 = 0;
			double acc5  = 0;		double acc15 = 0;
			double acc6  = 0;		double acc16 = 0;
			double acc7  = 0;		double acc17 = 0;
			double acc8  = 0;		double acc18 = 0;
			double acc9  = 0;		double acc19 = 0;
			double acc10 = 0;		double acc20 = 0;

			double in1  = 0;		double in11 = 0;
			double in2  = 0;		double in12 = 0;
			double in3  = 0;		double in13 = 0;
			double in4  = 0;		double in14 = 0;
			double in5  = 0;		double in15 = 0;
			double in6  = 0;		double in16 = 0;
			double in7  = 0;		double in17 = 0;
			double in8  = 0;		double in18 = 0;
			double in9  = 0;		double in19 = 0;
			double in10 = 0;		double in20 = 0;

			double sumNewton = a[0], 	sumDeriv = 0;
			double xpwrNewton = n, 		xpwrDeriv = 1;

			for (j = 1; j < n-9; j+=10)
			{
				in1 += xpwrDeriv *
					a[j] *(double)(j);
				in2 += xpwrDeriv *approx *
					a[j] * (double)(j+1);

				in3 += xpwrDeriv *approx *approx *
					a[j] * (double)(j+2);

				in4 += xpwrDeriv *approx *approx *approx *
					a[j] * (double)(j+3);

				in5 += xpwrDeriv *approx *approx *approx *approx *
					a[j] * (double)(j+4);

				in6 += xpwrDeriv *approx *approx *approx *approx *approx *
					a[j] * (double)(j+5);

				in7 += xpwrDeriv *approx *approx *approx *approx *approx *approx *
					a[j] * (double)(j+6);

				in8 += xpwrDeriv *approx *approx *approx *approx *approx *approx *approx *
					a[j] * (double)(j+7);

				in9 += xpwrDeriv *approx *approx *approx *approx *approx *approx *approx *approx *
					a[j] * (double)(j+8);

				in10 += xpwrDeriv *approx *approx *approx *approx *approx *approx *approx *approx *approx *
					a[j] * (double)(j+9);

				xpwrDeriv *= approx *approx *approx *approx *approx *approx *approx *approx *approx * approx;

				acc1 += in1;
				acc2 += in2;
				acc3 += in3;
				acc4 += in4;
				acc5 += in5;
				acc6 += in6;
				acc7 += in7;
				acc8 += in8;
				acc9 += in9;
				acc10 += in10;

				in11 += xpwrNewton * a[j];

				in12 += xpwrNewton*approx * a[j];

				in13 += xpwrNewton*approx *approx * a[j];

				in14 += xpwrNewton*approx *approx *approx * a[j];

				in15 += xpwrNewton*approx *approx *approx *approx * a[j];

				in16 += xpwrNewton*approx *approx *approx *approx *approx * a[j];

				in17 += xpwrNewton*approx *approx *approx *approx *approx *approx * a[j];

				in18 += xpwrNewton*approx *approx *approx *approx *approx *approx *approx * a[j];

				in19 += xpwrNewton*approx *approx *approx *approx *approx *approx *approx *approx * a[j];

				in20 += xpwrNewton*approx *approx *approx *approx *approx *approx *approx *approx *approx * a[j];

				xpwrNewton *=approx *approx *approx *approx *approx *approx *approx *approx *approx * approx;

				acc11 += in11;
				acc12 += in12;
				acc13 += in13;
				acc14 += in14;
				acc15 += in15;
				acc16 += in16;
				acc17 += in17;
				acc18 += in18;
				acc19 += in19;
				acc20 += in20;
				
			}

			//finish up
			for (; j < n+1; j++)
			{
				sumDeriv += a[j] * xpwrDeriv * (double)j;
				xpwrDeriv *= approx;

				sumNewton += a[j] * xpwrNewton;
				xpwrNewton *= approx;
			}

			//sum in dest
			check = sumNewton + acc1 + acc2 + acc3 + acc4 + acc5 + acc6 + acc7 + acc8 + acc9 + acc10;
			deriv = sumDeriv + acc11 + acc12 + acc13 + acc14 + acc15 + acc16 + acc17 + acc18 + acc19 +acc20;

			if (deriv == 0.0) check = n-1e-16;
			else check = approx + (double)(-1.0) * check/deriv;

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
