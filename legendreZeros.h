#ifndef NEW_LEGENDREZEROS_H
#define NEW_LEGENDREZEROS_H

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <iostream>


int getLegendreCoeff(double* coeff, int n);
int getLegendre(double* zero, double* weights, double* coeff, int n);
int errorCheck(double* coeff, int n);
void evaluate(double &approx, double &func, double &deriv,int n);
#endif