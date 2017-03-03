#include "stdafx.h"
#include "container.h"

counter::counter(bool f){
	cout<<"正在初始化四维数组"<<endl;
	flag = f;
	if (f)
	{	
		array.resize(20);
		prob.resize(20);
		aprob.resize(20);
		for (int i=0;i<20;i++)
		{
			array[i].resize(1296);
			
			cout<<"第一维初始化结束！"<<endl;
			for (int j=0;j<1296;j++)
			{
				array[i][j].resize(20);
				cout<<"第二维初始化结束！"<<endl;
				for (int k=0;k<20;k++){
					array[i][j][k].resize(1296);	
					cout<<"第三维初始化结束！"<<endl;
					for(int l=0;l<1296;l++)
						array[i][j][k][l] = 1;
				}
			}			
		}
	}else{
		array.resize(20);
		prob.resize(20);
		aprob.resize(20);
		for (int i=0;i<20;i++)
		{
			array[i].resize(15);
			prob[i].resize(15);
			aprob[i].resize(15);
			cout<<"第一维初始化结束！"<<endl;
			for (int j=0;j<15;j++)
			{
				array[i][j].resize(20);
				prob[i][j].resize(20);
				aprob[i][j].resize(20);
				cout<<"第二维初始化结束！"<<endl;
				for (int k=0;k<20;k++){
					array[i][j][k].resize(15);	
					prob[i][j][k].resize(15);
					aprob[i][j][k].resize(15);
					cout<<"第三维初始化结束！"<<endl;
					for(int l=0;l<15;l++)
					array[i][j][k][l] = 1;
				}
			}			
		}
	}
}

void counter::add(int x1, int y1, int x2, int y2){
	if(flag)
		++array[x1][y1][x2][y2]; 
	else
		++array[x1][y1-1][x2][y2-1];
}

void counter::comput(){
	cout<<"正在计算概率矩阵"<<endl;
	if (flag)//区别字母表的定义方案一、二
	{
		for (int i=0;i<20;i++)
			for (int j=0;j<1296;j++){
				double sum1 = 0.0;
				for (int k=0;k<20;k++)
				{
					double sum = 0.0;
					for(int l=0;l<1296;l++) sum = sum + array[i][j][k][l];
					for(int l=0;l<1296;l++) prob[i][j][k][l] = array[i][j][k][l] / sum;
					sum1 = sum1 + sum;
				}
				for(int k=0;k<20;k++)
					for(int l=0;l<1296;l++)
						aprob[i][j][k][l] = array[i][j][k][l] / sum1;
			}
	}
	else{
		for (int i=0;i<20;i++)
			for (int j=0;j<15;j++){
				double sum1 = 0.0;
				for (int k=0;k<20;k++)
				{
					double sum = 0.0;
					for(int l=0;l<15;l++) sum = sum + array[i][j][k][l];
					for(int l=0;l<15;l++) prob[i][j][k][l] = array[i][j][k][l] / sum;
					sum1 = sum1 + sum;
				}
				for(int k=0;k<20;k++)
					for(int l=0;l<15;l++)
						aprob[i][j][k][l] = array[i][j][k][l] / sum1;
			}
	}
}

void counter::print(ofstream & of){
	cout<<"正在输出第一个概率矩阵"<<endl;
	if(flag){
	for (int i=0;i<20;i++)
	{
		for (int j=0;j<1296;j++)
		{
			for (int k=0;k<20;k++)
			{
				for (int l=0;l<1296;l++) of<<prob[i][j][k][l]<<"  ";
				of<<endl;
			}
			of<<endl;
		}
		of<<endl;
	}
	}
	else{
		for (int i=0;i<20;i++)
		{
			for (int j=0;j<15;j++)
			{
				for (int k=0;k<20;k++)
				{
					for (int l=0;l<15;l++) of<<prob[i][j][k][l]<<"  ";
					of<<endl;
				}
				of<<endl;
			}
			of<<endl;
		}
	}
	cout<<"正在输出第二个概率矩阵"<<endl;
	if(flag){
		for (int i=0;i<20;i++)
		{
			for (int j=0;j<1296;j++)
			{
				for (int k=0;k<20;k++)
				{
					for (int l=0;l<1296;l++) of<<aprob[i][j][k][l]<<"  ";
					of<<endl;
				}
				of<<endl;
			}
			of<<endl;
		}
	}
	else{
		for (int i=0;i<20;i++)
		{
			for (int j=0;j<15;j++)
			{
				for (int k=0;k<20;k++)
				{
					for (int l=0;l<15;l++) of<<aprob[i][j][k][l]<<"  ";
					of<<endl;
				}
				of<<endl;
			}
			of<<endl;
		}
	}
}