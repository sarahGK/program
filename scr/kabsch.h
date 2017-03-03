#ifndef _KABSCH_H
#define _KABSCH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include "rd_pdb.h"
#include "eigen.h"
using namespace std;

#ifndef M_PI
#define M_PI          3.14159265358979323846
#endif
#define     EQN_EPS     1.e-9
#define	    IsZero(x)	((x) > -EQN_EPS && (x) < EQN_EPS)
#define cross(a,b,c) (a[0] = (b[1] * c[2]) - (b[2] * c[1]), \
		      a[1] = (b[2] * c[0]) - (b[0] * c[2]), \
                      a[2] = (b[0] * c[1]) - (b[1] * c[0]))
//代表三维空间中的一个点
//class Point
//{
//public:
//	Point();
//public:
//	double x,y,z;
//};
void matrix_transpose (double result[3][3], double a[3][3]);
double kabsch(const vector<Point>& strBuf1, const vector<Point>& strBuf2, double uu[3][3],double t[3]);
void matrix_multiply (double m[3][3], double a[3][3], double b[3][3]);
void eigen_values (double m[3][3], double values[3], double vectors[3][3]);
void my_eigen_values(double m[3][3], double values[3], double vectors[3][3]);//自己实现的球特征值和特征向量的方法
void transformpoint (double p_new[3], double m[3][3], double p[3]);
void transformpoint(Point& pnew,double m[3][3],Point p);
void normalise (double a[3]);
int cubic_roots (double c[4], double s[3]);
double cbrt(double q);//暂时认为该函数是计算立方根的
#endif