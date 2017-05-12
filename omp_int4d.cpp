// 4d integration
// g++ -o integrate test_integrate_feynman.c -lm -lrt
#define TIME_CODE 1
#define CPG 2.9
#define GIG 1000000000
#define PI 3.1415926535897932384626

#include "new_legendreZeros.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

double integrateGaussian4d(int n, double m);
void init_all(double* lZ, double* lWeights, double* lCoeff, double* sinVals, int n, double m);
void cleanUp(double* lZ, double* lWeights, double* lCoeff, double* sinVals, int n);
void sinCalc(double *sinVals, double *vals, int size);

int main(int argc, char** argv)
{

	int n = 100;
	double m = 0.01;

	if (argc > 1) n = atoi(argv[1]);

	struct timespec time1, time2, elapsed_integrate;
	struct timespec diff(struct timespec start, struct timespec end);

	printf("\nResult: %.15f\n", integrateGaussian4d(n,m));

	printf("\n");
	return 0;
}

// Gaussian
double integrateGaussian4d(int n, double m)
{
	// lZ = table of legendreZeros
	// lWeights = table of legendreWeights
	// lCoeff = table of legendre coefficients
	// sinVals = table used to store calculated value of zeros
	double* lZ;
	double* lWeights;
	double* lCoeff;
	double* sinVals;

	//printf("Done Declaring\n");
	int i, j, k, l;

#if TIME_CODE
	struct timespec time1, time2, elapsed;
	struct timespec diff(struct timespec start, struct timespec end);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif


	lZ = new double[n+1];
	lWeights = new double[n+1];	
	lCoeff = new double[n+1];
	sinVals = new double[n];

#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed = diff(time1, time2);
	double cputime1 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
	printf("\n--- Allocation time: %0.6f (msec)", cputime1/1000000);
#endif

#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	init_all(lZ, lWeights, lCoeff, sinVals, n, m);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed = diff(time1, time2);
	double cputime2 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
	printf("\n--- Initialization time: %0.6f (msec)", cputime2/1000000);
#endif


#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	double sum = 0;
	

/*********************************


			COMPUTATION


*********************************/

	double m2 = m*m;
#pragma omp parallel omp
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for(k = 0; k < n; k++)			
			{
				for (l = 0; l < n; l++)
				{
					sum += 1.0/pow((sinVals[i] + sinVals[j] + sinVals[k] + sinVals[l] + m2),2.0) * lWeights[i] * lWeights[j] * lWeights[k] * lWeights[l];
				}
			}
		}
	}
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed = diff(time1, time2);
	double cputime3 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
	printf("\n--- DotProduct time: %0.6f (msec)", cputime3/1000000);
#endif

#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	cleanUp(lZ, lCoeff, lWeights, sinVals, n);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed = diff(time1, time2);
	double cputime4 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
	printf("\nClean up time: %0.6f (msec)", cputime4/1000000);
#endif

	printf("\n=== Total time: %0.6f (msec)", (cputime1+cputime2+cputime3+cputime4)/1000000);

	return sum;
}

void init_all(double* lZ,
			 double* lCoeff,
			 double* lWeights,
			 double* sinVals,
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
	getLegendre(lZ, lWeights, lCoeff, n);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed_cpu = diff(time1, time2);
	double cputime = (double)((double)(CPG)*(double)(GIG * elapsed_cpu.tv_sec + elapsed_cpu.tv_nsec));
	printf("\n----- Legendre Init time: %0.6f (msec)", cputime/1000000);
#endif



// ---------------------------
// Initialize Legendre Weights
// ---------------------------
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	sinCalc(sinVals,lZ,n);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed_cpu = diff(time1, time2);
	cputime = (double)((double)(CPG)*(double)(GIG * elapsed_cpu.tv_sec + elapsed_cpu.tv_nsec));
	printf("\n----- Weight Init time: %0.6f (msec)", cputime/1000000);
#endif
}


void sinCalc(double *sinVals, double *vals, int size)
{
	for(int i = 0 ; i < size; i++)
	{
		sinVals[i] = sin(PI * vals[i]/2) * sin(PI * vals[i]/2);
	}
}

void cleanUp(double* lZ,
			 double* lCoeff,
			 double* lWeights,
			 double* sinVals,
			 int n)
{
	int i, j, k;
	delete[] lZ;
	delete[] lCoeff;
	delete[] lWeights;
	delete[] sinVals;
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