#define TIME_CODE 1
#define CPG 2.4
#define GIG 1000000000
#define PI 3.1415926535897932384626

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double trapezoidIntegration4D(int size, double m);
void sinCalc(double *pvals, double *vals, int size);

int main(int argc, char** argv)
{
	#if TIME_CODE
		struct timespec time1, time2, elapsed;
		struct timespec diff(struct timespec start, struct timespec end);
	#endif

	for(int i = 10; i <= 200; i += 10)
	{
		#if TIME_CODE
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
		#endif
		double result = trapezoidIntegration4D(i,0.01);
		#if TIME_CODE
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
			elapsed = diff(time1, time2);
			double cputime1 = (double)((double)(CPG)*(double)(GIG * elapsed.tv_sec + elapsed.tv_nsec));
			printf("\n--- Computation time for size %d: %0.6f (msec)\n",i, cputime1/1000000);
		#endif
		printf("\nResult for size : %0.6f\n", result);
	}
}

double trapezoidIntegration4D(int size, double m)
{
	double output = 0;
	double *values, *pvals;
	values = new double[size];
	pvals = new double[size];
	
	values[0] = -1;
	double difference = 2.0/size;
	for(int i = 0; i < size ; i++)
		values[i+1] = values[i]+ difference;


	sinCalc(pvals,values,size);

	double m2 = m*m;
	for(int i = 0; i < size; i++)
	{
		double mult = 1;
		if(i == 0 || i == size -1) mult /=2.0;
		for(int j = 0; j < size; j++)
		{
			if(j == 0 || j == size -1) mult /=2.0;
			for(int k = 0; k < size ; k++)
			{
				if(k == 0 || k == size -1) mult /=2.0;
				for(int l = 0; l < size ; l++)
				{
					if(l == 0 || l == size -1) mult /=2.0;
					output += 1/((pvals[i]+pvals[j]+pvals[k]+pvals[l]+m2) * (pvals[i]+pvals[j]+pvals[k]+pvals[l]+m2)) * mult;
				}
			}
				
		}
	}
	return output;

}

void sinCalc(double *pvals, double *vals, int size)
{
	for(int i = 0 ; i < size; i++)
	{
		pvals[i] = sin(PI * vals[i]/2) * sin(PI * vals[i]/2);
	}
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