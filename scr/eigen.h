/*************************************************

*************************************************/
#pragma once

#include <iostream>
#include "math.h"

using namespace std;

int eejcb(double * a,int n,double *v, double eps,int jt);
//void eastrq(double *a,int n,double *q,double *b,double *c);
void SortEigen(double mu[3],double a[3][3]);