/******************************************
 * ���������Ļ�������
 ******************************************/
#ifndef _VECTOR_H
#define _VECTOR_H
#include "rd_pdb.h"
#define EPS 1e-5
#define PIHALF          1.570796
#define PI              3.141593
#define TWOPI           6.283185
#define RADIAN          57.29578//һ��������Ӧ�Ķ���

double Distance(const Point p1,const Point p2); //����ռ���������֮��ľ���
Point Diff(const Point p1, const Point p2); //�����������Ӧ�Ķ�λʸ�������p2ָ��p1
Point Add(const Point p1,const Point p2);//�����ӷ�
double Dot(const Point p1, const Point p2);//����������λʸ���ĵ��
double Vector_Len(const Point p);			//���㶨λʸ���ĳ���
Point Vector_Stre(const Point p,double lam);//��λʸ��������
Point Cross(const Point p1, const Point p2);//����������λʸ�������
double Dihedralangle(const Point p1, const Point p2, const Point p3, const Point p4);//�����ĸ����γɵ������(-180,180)
double Atan2(double y, double x);//����������õķ�tan�ǣ����㷨
double Bond_Angle(const Point p1,const Point p2,const Point p3);//��������ԭ���γɵļ���
Point Vector_Rotate(const Point a,const Point z,const double theta);//��ʸ��a����ĳһ����ֱ��������z��ת�Ƚ�(�Ƕ�-180-180)֮��õ�����ʸ��
Point LastXYZ(const Point p1,const Point p2,const Point p3,const double theta,const double tao,const double bond_len);//��������ԭ�������������Լ�����(p2,p3,p4)�ͼ���(p3,p4)������ĸ��������p4
Point PointInLine(const Point p3,const Point z,const double bond_len);
#endif