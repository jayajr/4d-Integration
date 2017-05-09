// 4d integration
// g++ -o integrate test_integrate_feynman.c -lm -lrt
#define TIME_CODE 1
#define CPG 3.3
#define GIG 1000000000

#include "gaussianElimination.c"
#include "legendreZeros.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

double integrateGaussian4d(int n, int option);
void f1(double**** fxyzt, double* lZ, int n, double m);
void init_all(double* lZ, double**** fxyzt, double**** w, double** gQ_SoE, double* gQ_Soln, int n, double m);
void init_fxyzt(double*lZ, double**** fxyzt, int n, double m);
void init_w(double*lZ, double**** w, int n);
void init_gQ_SoE(double* lZ, double** gQ_SoE, int n);
void init_gQ_Soln(double* gQ_Soln, int n);
void cleanUp(double* lZ, double**** fxyzt, double**** w, double** gQ_SoE, double* gQ_Soln, int n);

int main(int argc, char** argv)
{

	int n = 30;
	double m = 0.01;

	struct timespec time1, time2, elapsed_integrate;
	struct timespec diff(struct timespec start, struct timespec end);

	printf("\nResult: %.15f\n", integrateGaussian4d(n,m));

	//for (n = 25; n < 26; n++)
	//{
		//for (m = 1/sqrt(10.0); m < 11; m += pow((double)(10.0), (double)(1.0/6.0)))
		//{
			//printf("%.15f, %.15f", m, integrateGaussian4d(n, m));
			//printf("\n");
		//}
	//}
	printf("\n");
	return 0;
}


// Integrates in 2 dimensions

// Input:
// n = order of legendre
// option = function to integrate

// Output:
// Numerical Integral of <option> function
// R = [-1, 1] x [-1, 1]

// Gaussian
double integrateGaussian4d(int n, int m)
{
	// lZ = table of legendreZeros
	// fxyzt = table of f(x, y)
	// w = table of gaussian weights
	double* lZ;
	double**** fxyzt;
	double**** w;
	double** gQ_SoE;
	double* gQ_Soln;

	//printf("Done Declaring\n");
	int i, j, k, l;

#if TIME_CODE
	struct timespec time1, time2, elapsed;
	struct timespec diff(struct timespec start, struct timespec end);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif


	lZ = new double[n+1];
	fxyzt = new double***[n+1];
		for(i = 0; i < n; i++)
		{
			fxyzt[i] = new double**[n];
			for (j = 0; j < n; j++)
			{
				fxyzt[i][j] = new double*[n];
				for (k = 0; k < n; k++)
				{
					fxyzt[i][j][k] = new double[n];
				}
			}
		}

	w = new double***[n+1];
		for(i = 0; i < n; i++)
		{
			w[i] = new double**[n];
			for (j = 0; j < n; j++)
			{
				w[i][j] = new double*[n];
				for (k = 0; k < n; k++)
				{
					w[i][j][k] = new double[n];
				}
			}
		}

	gQ_SoE = new double*[n+1];
		for(i = 0; i < n; i++) gQ_SoE[i] = new double[n];
	gQ_Soln = new double[n+1];

#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed = diff(time1, time2);
	double cputime1 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
	printf("\n--- Allocation time: %0.6f (msec)\n", cputime1/1000000);
#endif

#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	//printf("Done Allocating\n");
	init_all(lZ, fxyzt, w, gQ_SoE, gQ_Soln, n, m);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed = diff(time1, time2);
	double cputime2 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
	printf("\n--- Initialization time: %0.6f (msec)\n", cputime2/1000000);
#endif
	//for (i = 0; i < n; i++) printf("%.15f\n", lZ[i]);
	//if (i < n-1) printf("\t")
	//printf("\n");
	//printf("Completed init_all\n");


#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	double sum = 0;
	
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for(k = 0; k < n; k++)			
			{
				for (l = 0; l < n; l++)
				{
					sum += fxyzt[i][j][k][l] * w[i][j][k][l];
				}
			}
		}
	}
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed = diff(time1, time2);
	double cputime3 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
	printf("\n--- DotProduct time: %0.6f (msec)\n", cputime3/1000000);
#endif


	//printf("Attempting cleanUp\n");
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	cleanUp(lZ, fxyzt, w, gQ_SoE, gQ_Soln, n);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed = diff(time1, time2);
	double cputime4 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
	printf("\nClean up time: %0.6f (msec)\n", cputime4/1000000);
#endif

	printf("\n=== Total time: %0.6f (msec)\n", (cputime1+cputime2+cputime3+cputime4)/1000000);

	//printf("Completed cleanUp\n");
	return sum;
}

void f1(double**** fxyzt, double* lZ, int n, double m)
{
	int i, j, k, l;
	for (i = 0; i < n; i++)
	{
		double t0 = pow(sin(M_PI*lZ[i]/2.0), 2.0);
		for (j = 0; j < n; j++)
		{
			double t1 = pow(sin(M_PI*lZ[j]/2.0), 2.0);
			for (k = 0; k < n; k++)
			{
				double t2 = pow(sin(M_PI*lZ[k]/2.0), 2.0);
				for (l = 0; l < n; l++)
				{
					double t3 = pow(sin(M_PI*lZ[l]/2.0), 2.0);
					double m2 = pow(m, 2.0);
					fxyzt[i][j][k][l] = 1.0/pow(t0 + t1 + t2 + t3 + m2, 2.0);
				}
			}
		}
	}
	return;
}


void init_all(double* lZ,
			 double**** fxyzt,
			 double**** w,
			 double** gQ_SoE,
			 double* gQ_Soln,
			 int n,
			 double m)
{

// --------------
// Legendre Zeros
// --------------
#if TIME_CODE
	struct timespec time1, time2, elapsed_cpu;
	struct timespec diff(struct timespec start, struct timespec end);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	getZeroes(lZ, gQ_Soln, n);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed_cpu = diff(time1, time2);
	double cputime = (double)((double)(CPG)*(double)(GIG * elapsed_cpu.tv_sec + elapsed_cpu.tv_nsec));
	printf("\n----- Legendre Init time: %0.6f (msec)\n", cputime/1000000);
#endif


// ---------------------------------------
// Gaussian Quadrature System of Equations
// ---------------------------------------
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	init_gQ_SoE(lZ, gQ_SoE, n);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed_cpu = diff(time1, time2);
	cputime = (double)((double)(CPG)*(double)(GIG * elapsed_cpu.tv_sec + elapsed_cpu.tv_nsec));
	printf("\n----- GQ_SoE Init time: %0.6f (msec)\n", cputime/1000000);
#endif


// ----------------------------
// Gaussian Quadrature Solution
// ----------------------------
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	init_gQ_Soln(gQ_Soln, n);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed_cpu = diff(time1, time2);
	cputime = (double)((double)(CPG)*(double)(GIG * elapsed_cpu.tv_sec + elapsed_cpu.tv_nsec));
	printf("\n----- GQ_Soln Init time: %0.6f (msec)\n", cputime/1000000);
#endif




// --------------------
// Gaussian Elimination
// --------------------
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	gaussianElimination(gQ_SoE, gQ_Soln, n);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed_cpu = diff(time1, time2);
	cputime = (double)((double)(CPG)*(double)(GIG * elapsed_cpu.tv_sec + elapsed_cpu.tv_nsec));
	printf("\n----- Gaussian Elim time: %0.6f (msec)\n", cputime/1000000);
#endif


// -------------------
// Initialize Function
// -------------------
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	init_fxyzt(lZ, fxyzt, n, m);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed_cpu = diff(time1, time2);
	cputime = (double)((double)(CPG)*(double)(GIG * elapsed_cpu.tv_sec + elapsed_cpu.tv_nsec));
	printf("\n----- Function Init time: %0.6f (msec)\n", cputime/1000000);
#endif


// ---------------------------
// Initialize Legendre Weights
// ---------------------------
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	init_w(gQ_Soln, w, n);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed_cpu = diff(time1, time2);
	cputime = (double)((double)(CPG)*(double)(GIG * elapsed_cpu.tv_sec + elapsed_cpu.tv_nsec));
	printf("\n----- Weight Init time: %0.6f (msec)\n", cputime/1000000);
#endif
}


void init_fxyzt(double* lZ, double**** fxyzt, int n, double m)
{
	f1(fxyzt, lZ, n, m);
}

void init_w(double* gQ_Soln, double**** w, int n)
{
	//	printf("gQSolni, gQSolnj, w:\n");
	int i, j, k, l;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < n; k++)
			{
				for (l = 0; l < n; l++)
				{
					w[i][j][k][l] = gQ_Soln[i] * gQ_Soln[j] * gQ_Soln[k] * gQ_Soln[l];
					//printf("%.15f, %.15f, %.15f\n", gQ_Soln[i], gQ_Soln[j], w[i][j]);
				}
			}
		}
	}
	return;
}

void init_gQ_SoE(double* lZ, double** gQ_SoE, int n)
{
	int i, j;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			gQ_SoE[i][j] = pow(lZ[j], (double)i);
		}
	}
	return;
}

void init_gQ_Soln(double* gQ_Soln, int n)
{
	int i;
	for(i = 0; i < n; i++)
	{
		if (i%2 == 0)
			gQ_Soln[i] = 2.0/(i+1);
		else
			gQ_Soln[i] = 0;
	}

	return;
}

void cleanUp(double* lZ,
			 double**** fxyzt,
			 double**** w,
			 double** gQ_SoE,
			 double* gQ_Soln,
			 int n)
{
	int i, j, k;
	delete[] lZ;
	
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < n; k++)
			{
				delete fxyzt[i][j][k];
			}
			delete fxyzt[i][j];
		}
		delete fxyzt[i];
	}
	delete[] fxyzt;
	
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < n; k++)
			{
				delete w[i][j][k];
			}
			delete w[i][j];
		}
		delete w[i];
	}
	delete[] w;

	for (i = 0; i < n; i++) delete[] gQ_SoE[i];
	delete[] gQ_SoE;

	delete[] gQ_Soln;
	return;
}

struct timespec diff(struct timespec start, struct timespec end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}