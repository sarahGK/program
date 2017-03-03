/*******************************************************
 * 在CE源程序上的KABSCH算法基础上所作的改动
 * 针对原方法中对于规则形状或者完全一样的形状
 * 不能很好处理的问题，使用了自定义的球特征值
 * 和特征向量的函数，但是根据基于gsl库的kabsch算法
 * 求得的结果需要排序
 ******************************************************/
#include "stdafx.h"
#include "kabsch.h"
#include "eigen.h"
#include "time.h"

//Point::Point()
//{
//	x=0;y=0;z=0;
//}
/*****************************************************
 * 两个点集在空间中最优重叠变换的KABSCH算法
 * 返回两个点集最优变换后的RMSD
 * 如果两个点集个数不一样，返回负值
 * 注意
 * 1. strBuf2中的点和旋转矩阵以及平移向量作用
 * 后，得到strBuf1中的点
 * 2. 每一个点的三维坐标需要看作行向量和旋转
 * 矩阵作用
 ****************************************************/
double kabsch(const vector<Point>& strBuf1, const vector<Point>& strBuf2, double uu[3][3],double t[3])
{

/*
 * Fit pairs of points using the method described in:
 *
 * `A discussion of the solution for the best rotation to relate
 *  two sets of vectors', W.Kabsch, Acta Cryst. (1978), A34, 827-828.
 *
 * The method generates the 4x4 matrix rot which will map the
 * coordinates y on to the coordinates x.
 *
 * The matrix is arranged such that rot[3][0], rot[3][1], rot[3][2]
 * are the translational components.
 *
 * The function returns the root mean square distance between
 * x and y.  The coordinates in x and y are not modified.
 */


  double rot[4][4];
  double xc[3], yc[3];
  double r[3][3], rtrans[3][3], rr[3][3];
  double mu[3], a[3][3], atrans[3][3], b[3][3], u[3][3];
  double yy[3];
  double d[20];
  double rmsd;
  int i, j, k;
  
  int n;

  for(i = 0; i < 3; i++) {
    xc[i] = 0.0;
    yc[i] = 0.0;
  } 
if(strBuf1.size()!=strBuf2.size())
return -1;
  n =(int)strBuf1.size();
  

  /* Find center of each set of coordinates. */

  for (i = 0; i < n; i++)
    {
      xc[0] += strBuf1[i].x;
      yc[0] += strBuf2[i].x;
      xc[1] += strBuf1[i].y;
      yc[1] += strBuf2[i].y;
      xc[2] += strBuf1[i].z;
      yc[2] += strBuf2[i].z;
      
      //for (j = 0; j < 3; j++)
      //{
      //  xc[j] += MX[i][j];
      //  yc[j] += MY[i][j];
      //}
    }

  for (j = 0; j < 3; j++)
    {
      xc[j] /= (double) n;
      yc[j] /= (double) n;
    }

  /*
   * Initialise and then fill the r matrix.
   * Note that centre is subtracted at this stage.
   */

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      r[i][j] = 0.0;

  for (k = n-1; k >=0; k--)
    {

      r[0][0] += (strBuf2[k].x - yc[0]) * (strBuf1[k].x - xc[0]);
      r[0][1] += (strBuf2[k].x - yc[0]) * (strBuf1[k].y - xc[1]);
      r[0][2] += (strBuf2[k].x - yc[0]) * (strBuf1[k].z - xc[2]);
      r[1][0] += (strBuf2[k].y - yc[1]) * (strBuf1[k].x - xc[0]);
      r[1][1] += (strBuf2[k].y - yc[1]) * (strBuf1[k].y - xc[1]);
      r[1][2] += (strBuf2[k].y - yc[1]) * (strBuf1[k].z - xc[2]);
      r[2][0] += (strBuf2[k].z - yc[2]) * (strBuf1[k].x - xc[0]);
      r[2][1] += (strBuf2[k].z - yc[2]) * (strBuf1[k].y - xc[1]);
      r[2][2] += (strBuf2[k].z - yc[2]) * (strBuf1[k].z - xc[2]);

      //for (i = 0; i < 3; i++)
      //{
      //  for (j = 0; j < 3; j++)
      //    {
      //      r[i][j] += (MY[k][i] - yc[i]) * (MX[k][j] - xc[j]);
      //    }
      //}

    }

  /* Generate the transpose of r and form rtrans x r */

  matrix_transpose (rtrans, r);

  matrix_multiply (rr, rtrans, r);


  ////输出矩阵rr
  //for(i=0;i<3;i++)
  //{
	 // for(j=0;j<3;j++)
		//  cout<<rr[i][j]<<'\t';
	 // cout<<endl;
  //}
  /*
   * Get the eigenvalues and vectors.
   * Reform a[2] as cross product of a[0] and a[1] to ensure
   * right handed system.
   */

  //增加
  for(i=0;i<3;i++)
  {
	  mu[i]=0;
	  for(j=0;j<3;j++)
		  a[i][j]=0;
  }
  //eigen_values (rr, mu, a);
  

  my_eigen_values(rr,mu,a);//使用自定义的求特征值和特征向量函数 特征值和特征向量的排列和原始方法不一致，所以需要排序  
  SortEigen(mu,a);

  cross (a[2], a[0], a[1]);

  /* Transform first two eigenvectors and normalise them. */

  for (i = 0; i < 2; i++)
    {
      transformpoint (b[i], r, a[i]);

      normalise (b[i]);
    }

  /* Make right handed set. */

  cross (b[2], b[0], b[1]);

  /* Form the rotation matrix. */

  matrix_transpose (atrans, a);

  matrix_multiply (u, b, atrans);

  /* Make rot the identity matrix. */

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      rot[i][j] = (i == j) ? 1. : 0.;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      rot[i][j] = u[i][j];

  /* Transform offset of y coordinates by the rotation. */

  transformpoint (yy, u, yc);

  /* Build translational component of rot from offsets. */

  for (i = 0; i < 3; i++)
    rot[3][i] = -yy[i] + xc[i];

  /* Figure out the rms deviation of the fitted coordinates. */

  d[0] = rot[0][0]; d[1] = rot[1][0]; d[2] = rot[2][0]; 
  d[3] = rot[0][1]; d[4] = rot[1][1]; d[5] = rot[2][1]; 
  d[6] = rot[0][2]; d[7] = rot[1][2]; d[8] = rot[2][2]; 
  d[9] = rot[3][0]; d[10] = rot[3][1]; d[11] = rot[3][2]; 

  for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
		  uu[i][j]=u[i][j];
  for(i=0;i<3;i++)
	  t[i]=rot[3][i];


  rmsd = 0.0;

  for (i = 0; i < n; i++)
    {
      double xt[3], yt[3];

      /* Translate the pairs of coordinates to origin. */

      xt[0] = strBuf1[i].x - xc[0];
      xt[1] = strBuf1[i].y - xc[1];
      xt[2] = strBuf1[i].z - xc[2];

      yt[0] = strBuf2[i].x - yc[0];
      yt[1] = strBuf2[i].y - yc[1];
      yt[2] = strBuf2[i].z - yc[2];

	  //for (j = 0; j < 3; j++)
      //{
      //  xt[j] = MX[i][j] - xc[j];
      //  yt[j] = MY[i][j] - yc[j];
      //}

      /* Transform the y coordinate by u. */

      transformpoint (yy, u, yt);
	  

      /* Accumulate the difference. */

      rmsd +=
	(yy[0] - xt[0]) * (yy[0] - xt[0]) +
	(yy[1] - xt[1]) * (yy[1] - xt[1]) +
	(yy[2] - xt[2]) * (yy[2] - xt[2]);
    }

  rmsd /= (double) n;

  rmsd = sqrt (rmsd);

  return rmsd;
}
void matrix_transpose (double result[3][3], double a[3][3])
{
  int i, j;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      result[i][j] = a[j][i];
}
void matrix_multiply (double m[3][3], double a[3][3], double b[3][3])
{
  int row, column, i;

  for (row = 0; row < 3; row++)
    for (column = 0; column < 3; column++)
      {
	m[row][column] = 0.0;
	for (i = 0; i < 3; i++)
	  m[row][column] += b[row][i] * a[i][column];
      }
}
void eigen_values (double m[3][3], double values[3], double vectors[3][3])
/*
 * Find the eigen values and vectors for the matrix m.
 */
{
  double a1 = m[0][0], b1 = m[0][1], c1 = m[0][2];
  double a2 = m[1][0], b2 = m[1][1], c2 = m[1][2];
  double a3 = m[2][0], b3 = m[2][1], c3 = m[2][2];
  double c[4];
  double x, y, z, norm;
  double roots[3], l;
  int nroots, iroot;
  int i;

  /*
   * Expanding the characteristic equation of a 3x3 matrix...
   * Maple tells us the following is the cubic in l
   */

  c[0] = a1 * (b2 * c3 - b3 * c2) + c1 * (a2 * b3 - a3 * b2) - b1 * (a2 * c3 - a3 * c2);
  c[1] = (a1 * (-b2 - c3) - b2 * c3 + b3 * c2 + c1 * a3 + b1 * a2);
  c[2] = (a1 + b2 + c3);
  c[3] = -1.;

  nroots = cubic_roots (c, roots);

  /* Degenerate roots are not returned individually. */

  for (i = 0; i < nroots; i++)
    {
      iroot = i > nroots ? nroots : i;

      l = roots[iroot];

      values[i] = l;

      /*
       * Find the eigen vectors by solving pairs of the
       * three simultaneous equations.  From `Mathematical Methods
       * in Science and Engineering', Heiding, pg.19
       *
       * Sometimes we get x = y = z = 0.0, so try the other two
       * pairs of equations and hope that one of them gives a solution.
       */

      x = b1 * c2 - (b2 - l) * c1;
      y = -((a1 - l) * c2 - a2 * c1);
      z = ((a1 - l) * (b2 - l) - a2 * b1);

      if (IsZero (x) && IsZero (y) && IsZero (z))
	{
	  x = b1 * (c3 - l) - b3 * c1;
	  y = -((a1 - l) * (c3 - l) - a3 * c1);
	  z = ((a1 - l) * b3 - a3 * b1);

	  if (IsZero (x) && IsZero (y) && IsZero (z))
	    {
	      x = (b2 - l) * (c3 - l) - b3 * c2;
	      y = -(a2 * (c3 - l) - a3 * c2);
	      z = (a2 * b3 - a3 * (b2 - l));
	      if (IsZero (x) && IsZero (y) && IsZero (z))
		{
		  printf ("eigen: no solution for eigen vector %d\n", i);
		  exit(0);			
		}
	    }
	}

      norm = sqrt (x * x + y * y + z * z);

      if (!IsZero (norm))
	{
	  vectors[i][0] = x / norm;
	  vectors[i][1] = y / norm;
	  vectors[i][2] = z / norm;
	}
    }
}
/*****************************************
 * 单个点的旋转
 * 注意点看作是行向量，经过矩阵和
 * m相乘后得到p_new
 * 如果把点看作是列向量需要把m转置
 *****************************************/
void transformpoint (double p_new[3], double m[3][3], double p[3])
{
  int i, column;

  for (i = 0; i < 3; i++)
    {
      p_new[i] = 0.0;
      for (column = 0; column < 3; column++)
	p_new[i] += m[column][i] * p[column];
    }
}
void transformpoint(Point& pnew,double m[3][3],Point p)
{
	double pp[3],p_new[3];
	pp[0]=p.x;
	pp[1]=p.y;
	pp[2]=p.z;
	transformpoint(p_new,m,pp);
	pnew.x=p_new[0];
	pnew.y=p_new[1];
	pnew.z=p_new[2];
}
void normalise (double a[3])
{
  double norm = sqrt (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

  if (norm > 1.e-9)
    {
      a[0] /= norm;
      a[1] /= norm;
      a[2] /= norm;
    }
}
int cubic_roots (double c[4], double s[3])
{
  int i, num;
  double sub;
  double A, B, C;
  double A2, p, q;
  double p3, D;

  /* normal form: x^3 + Ax^2 + Bx + C = 0 */

  A = c[2] / c[3];
  B = c[1] / c[3];
  C = c[0] / c[3];

  /*  substitute x = y - A/3 to eliminate quadric term:
     x^3 +px + q = 0 */

  A2 = A * A;
  p = 1.0 / 3 * (-1.0 / 3 * A2 + B);
  q = 1.0 / 2 * (2.0 / 27 * A * A2 - 1.0 / 3 * A * B + C);

  /* use Cardano's formula */

  p3 = p * p * p;
  D = q * q + p3;

  if (IsZero (D))
    {				/* one triple solution */
      if (IsZero (q))
	{
	  s[0] = 0;
	  num = 1;
	}
      else
	{			/* one single and one double solution */
	  double u = cbrt (-q);
	  s[0] = 2 * u;
	  s[1] = -u;
	  num = 2;
	}
    }
  else if (D < 0)
    {				/* Casus irreducibilis: three real solutions */

      double phi = 1.0 / 3 * acos (-q / sqrt (-p3));
      double t = 2 * sqrt (-p);

      s[0] = t * cos (phi);
      s[1] = -t * cos (phi + M_PI / 3);
      s[2] = -t * cos (phi - M_PI / 3);
      num = 3;
    }
  else
    {				/* one real solution */

      double sqrt_D = sqrt (D);
      double u = cbrt (sqrt_D - q);
	  double temp=cbrt(sqrt_D+q);
      double v = -cbrt (sqrt_D + q);

      s[0] = u + v;
      num = 1;
    }

  /* resubstitute */

  sub = 1.0 / 3 * A;

  for (i = 0; i < num; ++i)
    s[i] -= sub;

  return num;
}
double cbrt(double q)
{
	double tempq;
	if(q<0)
		tempq=-1.0*q;
	else
		tempq=q;
	double re=pow(tempq,1.0/3);
	if(q<0)
		re=-1.0*re;
	return re;
}
void my_eigen_values(double m[3][3], double values[3], double vectors[3][3])
{
	double mm[3][3],vt[3][3],norm;
	int i,j;
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			mm[i][j]=m[i][j];
	int result=eejcb(&mm[0][0],3,&vt[0][0],1e-5,100);
	if(result<0)
	{
		cout<<"无法求特征值"<<endl;
		exit(0);
	}
	for(i=0;i<3;i++)
		values[i]=mm[i][i];
	for(i=0;i<3;i++)
	{
		norm=0;
		for(j=0;j<3;j++)
			norm+=vt[j][i]*vt[j][i];
		norm=sqrt(norm);
		for(j=0;j<3;j++)
			vt[j][i]/=norm;
	}
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			vectors[i][j]=vt[j][i];
}

