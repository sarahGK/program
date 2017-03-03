#include "stdafx.h"
#include "vector.h"


double Distance(const Point p1,const Point p2)
{
	double x2,y2,z2;
	x2=(p1.x-p2.x)*(p1.x-p2.x);
	y2=(p1.y-p2.y)*(p1.y-p2.y);
	z2=(p1.z-p2.z)*(p1.z-p2.z);
	return sqrt(x2+y2+z2);
}
Point Diff(const Point p1, const Point p2)
{
	Point result;
	result.x=p1.x-p2.x;
	result.y=p1.y-p2.y;
	result.z=p1.z-p2.z;
	return result;
}
Point Add(const Point p1,const Point p2)
{
	Point result;
	result.x=p1.x+p2.x;
	result.y=p1.y+p2.y;
	result.z=p1.z+p2.z;
	return result;
}
double Dot(const Point p1, const Point p2)
{
	double result;
	result=p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
	return result;
}
Point Cross(const Point p1, const Point p2)
{
	Point result;
	result.x=p1.y*p2.z-p2.y*p1.z;
	result.y=p1.z*p2.x-p2.z*p1.x;
	result.z=p1.x*p2.y-p2.x*p1.y;
	return result;
}
/*******************************************************
 * 计算四个点形成的两面角
 * 说明：两面角的取值为(-180,180)之间
 * 计算它的关键问题是确定该两面角是正的还是负的
 * 以及他是锐角还是钝角
 * 以下方法中atan2(v,u)中的v的正负决定角度的正负
 * u的正负决定该角是锐角还是钝角
 ********************************************************/
double Dihedralangle(const Point v1, const Point v2, const Point v3, const Point v4)
{
	double Result, u, v;
  Point v12, v43, x, y, z, p;

  v12=Diff(v1, v2);
  v43=Diff(v4, v3);
  z=Diff(v2, v3);
  p=Cross(z, v12);
  x=Cross(z, v43);
  y=Cross(z, x);
  u = Dot(x, x);
  v = Dot(y, y);
  Result = 360.0;
  if (u <= 0.0 || v <= 0.0)
    return Result;
  u = Dot(p, x) / sqrt(u);
  v = Dot(p, y) / sqrt(v);
  if (u != 0.0 || v != 0.0)
    return (Atan2(v, u) * RADIAN);
  return Result;
}
double Atan2(double y, double x)
{
  double z;

  if (x != 0.0)
    z = atan(y / x);
  else if (y > 0.0)
    z = PIHALF;
  else if (y < 0.0)
    z = -PIHALF;
  else
    z = TWOPI;
  if (x >= 0.0)
    return z;
  if (y > 0.0)
    z += PI;
  else
    z -= PI;
  return z;
}
double Bond_Angle(const Point p1,const Point p2,const Point p3)
{
	Point a,b;
	a=Diff(p1,p2);
	b=Diff(p3,p2);
	double c=Dot(a,b);
	double result,len;
	len=sqrt(Dot(a,a))*sqrt(Dot(b,b));
	result=acos(c/len);
	return result*RADIAN;
}
double Vector_Len(const Point p)
{
	return sqrt(Dot(p,p));
}
Point Vector_Stre(const Point p,double lam)
{
	Point a;
	a.x=p.x*lam;
	a.y=p.y*lam;
	a.z=p.z*lam;
	return a;
}
Point Vector_Rotate(const Point a,const Point z,const double theta)
{
	Point c,cc,b;
	double lamba,theta1;
	if(fabs(theta-90)<EPS)
		b=Cross(z,a);
	else if(fabs(theta+90)<EPS)
		b=Cross(a,z);
	else if((theta>=0)&&(theta<90))
	{
		c=Cross(z,a);
		lamba=Vector_Len(a)*tan(theta/RADIAN)/Vector_Len(c);
		cc=Vector_Stre(c,lamba);
		b=Add(a,cc);
	}
	else if((theta>90)&&(theta<=180))
	{
		c=Cross(z,a);
		lamba=Vector_Len(a)*tan((180-theta)/RADIAN)/Vector_Len(c);
		cc=Vector_Stre(c,lamba);
		b=Diff(cc,a);
	}
	else if((theta<0)&&(theta>-90))
	{
		c=Cross(a,z);
		lamba=Vector_Len(a)*tan(-1.0*theta/RADIAN)/Vector_Len(c);
		cc=Vector_Stre(c,lamba);
		b=Add(a,cc);
	}
	else if((theta<-90)&&(theta>=-180))
	{
		c=Cross(a,z);
		lamba=Vector_Len(a)*tan((180+theta)/RADIAN)/Vector_Len(c);
		cc=Vector_Stre(c,lamba);
		b=Diff(cc,a);
	}
	else
	{
		if(theta>180)
		{
			theta1=theta-360;
			b=Vector_Rotate(a,z,theta1);
		}
		else
		{
			theta1=theta+360;
			b=Vector_Rotate(a,z,theta1);
		}
	}
	return b;
}
Point LastXYZ(const Point p1,const Point p2,const Point p3,const double theta,const double tao,const double bond_len)
{
	Point aa,a,z,result;
	Point p21,p23;
	p21=Diff(p1,p2);
	p23=Diff(p3,p2);
	aa=Cross(p23,p21);
	a=Vector_Rotate(aa,p23,theta);
	z=Vector_Rotate(p23,a,180-tao);
	result=PointInLine(p3,z,bond_len);
	return result;
}
/***********************************************************
 * 计算直线(经过点p3,方向矢量为z)上与p3距离为bond_len
 * 且与z同向的点
 **********************************************************/
Point PointInLine(const Point p3,const Point z,const double bond_len)
{
	double t;
	Point a,b;
	t=bond_len/Vector_Len(z);
	b=Vector_Stre(z,t);
	a=Add(p3,b);
	return a;
}


