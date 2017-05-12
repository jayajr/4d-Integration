// 4d integration
// g++ -o integrate test_integrate_feynman.c -lm -lrt
#define TIME_CODE 1
#define CPG 1.6
#define GIG 1000000000
#define PI 3.1415926535897932384626
#define M2 0.0001

#include "gpu_legendreZeros.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/execution_policy.h>


double integrateGaussian4d(int n, double m);
void init_all(double* lZ, double* lWeights, double* lCoeff, double* sinVals, int n, double m);
void cleanUp(double* lZ, double* lWeights, double* lCoeff, double* sinVals, int n);
void sinCalc(double *sinVals, double *vals, int size);

struct permuteFunctor
{
	const double* sinVals_f;
	const double* weights_f;
	double* sumVec_f;
	int n;

	permuteFunctor(thrust::device_vector<double> const& sinV,
				   thrust::device_vector<double> const& gQV,
				   thrust::device_vector<double> & sumV)
	{
		sinVals_f = thrust::raw_pointer_cast(sinV.data());
		weights_f = thrust::raw_pointer_cast(gQV.data());
		sumVec_f = thrust::raw_pointer_cast(sumV.data());
		n = sinV.size();
	}



	__device__
	void operator()(int x)
	{
		int nn = (n*n);
		int nnn = (n*n*n);

		int l = (x % (n));
		int k = (x % (nn)  / (n));
		int j = (x % (nnn) / (nn));
		int i = (x         / (nnn));

		double sumSin = sinVals_f[i] + sinVals_f[j] + sinVals_f[k] + sinVals_f[l] + M2;
		double prodGQ = weights_f[i] * weights_f[j] * weights_f[k] * weights_f[l];

		sumVec_f[(i*(nnn)) + (j*(nn)) + (k*(n)) + l] = (1.0/(sumSin*sumSin)*prodGQ);

	}
};

int main(int argc, char** argv)
{

	int n = 100;
	double m = 0.01;

	if (argc > 1) n = atoi(argv[1]);

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
	printf("\n--- Allocation time: %0.6f (msec)\n", cputime1/1000000);
#endif

#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	init_all(lZ, lWeights, lCoeff, sinVals, n, m);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed = diff(time1, time2);
	double cputime2 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
	printf("\n--- Initialization time: %0.6f (msec)\n", cputime2/1000000);
#endif



#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	double sum = 0;
	

/*********************************


				GPU


*********************************/


	thrust::device_vector<double> d_sinVals(n);
	thrust::device_vector<double> d_Weights(n);
	thrust::device_vector<double> d_sumVals(n*n*n*n, 0);

	for (int i = 0; i < n; i ++)
	{
		d_sinVals[i] = sinVals[i];
		d_Weights[i] = lWeights[i];
	}

	thrust::for_each_n(
		thrust::device,
		thrust::counting_iterator<int>(0),
		(n*n*n*n),
		permuteFunctor(d_sinVals, d_Weights, d_sumVals));


	sum = thrust::reduce(d_sumVals.begin(), d_sumVals.end(), (int) 0, thrust::plus<int>());


#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed = diff(time1, time2);
	double cputime3 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
	printf("\n--- DotProduct time: %0.6f (msec)\n", cputime3/1000000);
#endif


#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
#endif
	cleanUp(lZ, lCoeff, lWeights, sinVals, n);
#if TIME_CODE
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	elapsed = diff(time1, time2);
	double cputime4 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
	printf("\nClean up time: %0.6f (msec)\n", cputime4/1000000);
#endif

	printf("\n=== Total time: %0.6f (msec)\n", (cputime1+cputime2+cputime3+cputime4)/1000000);

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
	printf("\n----- Legendre Init time: %0.6f (msec)\n", cputime/1000000);
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
	printf("\n----- Weight Init time: %0.6f (msec)\n", cputime/1000000);
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
			 double* lWeights,
			 double* lCoeff,
			 double* sinVals,
			 int n)
{
	int i;
	delete[] lZ;
	delete[] lWeights;
	delete[] lCoeff;
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