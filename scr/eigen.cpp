#include "stdafx.h"
#include "eigen.h"

/****************************************************
 * 求实对称矩阵的特征值和特征向量的雅可比法
 * 参数:
 * a: 数组,体积为n*n,存放n阶实对称矩阵,返回时对
 * 角线上存放n个特征值
 * n:矩阵的阶数
 * v:数组,体积为n*n,存放特征向量,其中第i列为与
 * lambai对应的特征向量
 * eps控制精度
 * jt最大的叠带次数
 * 返回值:小于0,表明叠带了jt次还没有找到特征值
 * 大于0,则正常
 * 参考书<<C常用算法程序集>>P205
 *****************************************************/
int eejcb(double * a,int n,double *v, double eps,int jt)
{
	int i,j,p,q,u,w,t,s,l;
	double fm,cn,sn,omega,y,x,d;
	l=1;

	for(i=0;i<=n-1;i++)
	{
		v[i*n+i]=1.0;
		for(j=0;j<=n-1;j++)
			if(i!=j)
				v[i*n+j]=0;
	}
	//while(l==1)
	while(l<jt)
	{
		fm=0;
		for(i=0;i<=n-1;i++)
			for(j=0;j<=n-1;j++)
			{
				d=fabs(a[i*n+j]);
				if((i!=j)&&(d>fm))
				{
					fm=d;
					p=i;
					q=j;
				}
			}
		if(fm<eps)
			return(1);
		if(l>jt)
			return(-1);
		l=l+1;
		u=p*n+q;w=p*n+p;t=q*n+p;s=q*n+q;
		x=-a[u];y=(a[s]-a[w])/2;
		omega=x/sqrt(x*x+y*y);
		if(y<0)
			omega=-omega;
		sn=1.0+sqrt(1.0-omega*omega);
		sn=omega/sqrt(2.0*sn);
		cn=sqrt(1.0-sn*sn);
		fm=a[w];
		a[w]=fm*cn*cn+a[s]*sn*sn+a[u]*omega;
		a[s]=fm*sn*sn+a[s]*cn*cn-a[u]*omega;
		a[u]=0;a[t]=0;
		for(j=0;j<=n-1;j++)
			if((j!=p)&&(j!=q))
			{
				u=p*n+j;
				w=q*n+j;
				fm=a[u];
				a[u]=fm*cn+a[w]*sn;
				a[w]=-fm*sn+a[w]*cn;
			}
		for(i=0;i<=n-1;i++)
			if((i!=p)&&(i!=q))
			{
				u=i*n+p;w=i*n+q;
				fm=a[u];
				a[u]=fm*cn+a[w]*sn;
				a[w]=-fm*sn+a[w]*cn;
			}
		for(i=0;i<=n-1;i++)
		{
			u=i*n+p;w=i*n+q;
			fm=v[u];
			v[u]=fm*cn+v[w]*sn;
			v[w]=-fm*sn+v[w]*cn;
		}
	}
	return(1);
}

void SortEigen(double mu[3],double a[3][3])
{
	int i,j,k;
	double temp,tempmu[3],tempa[3][3];
	double ej,ek;

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
			tempa[i][j]=a[i][j];//根据gsl的kabsch算法，取转置?
		tempmu[i]=mu[i];
	}
	for(i=0;i<2;i++)
	{
		k=i;
		ek=tempmu[i];
		for(j=i+1;j<3;j++)
		{
			ej=tempmu[j];
			if(ej>ek)
			{
				k=j;
				ek=ej;
			}
		}
		if(k!=i)
		{
			temp=tempmu[i];tempmu[i]=tempmu[k];tempmu[k]=temp;
			//交换i,k列
			/*temp=tempa[0][i];tempa[0][i]=tempa[0][k];tempa[0][k]=temp;
			temp=tempa[1][i];tempa[1][i]=tempa[1][k];tempa[1][k]=temp;
			temp=tempa[2][i];tempa[2][i]=tempa[2][k];tempa[2][k]=temp;*/
			//交换i,k行
			temp=tempa[i][0];tempa[i][0]=tempa[k][0];tempa[k][0]=temp;
			temp=tempa[i][1];tempa[i][1]=tempa[k][1];tempa[k][1]=temp;
			temp=tempa[i][2];tempa[i][2]=tempa[k][2];tempa[k][2]=temp;
		}
	}

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
			a[i][j]=tempa[i][j];
		mu[i]=tempmu[i];
	}
}
