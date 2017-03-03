/******************************************
 * 定义向量的基本运算
 ******************************************/
#ifndef _VECTOR_H
#define _VECTOR_H
#include "rd_pdb.h"
#define EPS 1e-5
#define PIHALF          1.570796
#define PI              3.141593
#define TWOPI           6.283185
#define RADIAN          57.29578//一弧度所对应的度数

double Distance(const Point p1,const Point p2); //计算空间中两个点之间的距离
Point Diff(const Point p1, const Point p2); //计算两个点对应的定位矢量方向从p2指向p1
Point Add(const Point p1,const Point p2);//向量加法
double Dot(const Point p1, const Point p2);//计算两个定位矢量的点积
double Vector_Len(const Point p);			//计算定位矢量的长度
Point Vector_Stre(const Point p,double lam);//定位矢量的数乘
Point Cross(const Point p1, const Point p2);//计算两个定位矢量的外积
double Dihedralangle(const Point p1, const Point p2, const Point p3, const Point p4);//计算四个点形成的两面角(-180,180)
double Atan2(double y, double x);//计算两面角用的反tan角，计算法
double Bond_Angle(const Point p1,const Point p2,const Point p3);//计算三个原子形成的键角
Point Vector_Rotate(const Point a,const Point z,const double theta);//求矢量a绕着某一个垂直于它的轴z旋转θ角(角度-180-180)之后得到的新矢量
Point LastXYZ(const Point p1,const Point p2,const Point p3,const double theta,const double tao,const double bond_len);//给定三个原子坐标和两面角以及键角(p2,p3,p4)和键长(p3,p4)，求第四个点的坐标p4
Point PointInLine(const Point p3,const Point z,const double bond_len);
#endif