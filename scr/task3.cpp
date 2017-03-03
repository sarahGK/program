// task3.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "globfuc.h"
#include "container.h"
#include "hr_Info.h"
#include "readin_Array.h"


void  main()
{
	ifstream datain;
	datain.open("E:\\talk\\biology\\task2\\alphabet.txt",ios::in);//打开hr的字母表文件
	if(!datain.is_open()){
		cout<<"can't open the file: E:\\talk\\biology\\task2\\alphabet.txt"<<endl;
	}else{
		int t1 = 0;
		string str1;	//记录前氨基酸的信息
		counter c2(false);//定义第二种字母表定义方案的计数器
		cout<<"开始读取文件的数据"<<endl;
		getline(datain,str1);//读取文件的第一行
		while(!datain.eof()){
			t1 ++;
			if(str1.length()>10){//如果读取的不是pdb文件名
				char n1 = str1.substr(0,1).c_str()[0];//氨基酸的字母表示
				int x1 = atoi(str1.substr(3,5).c_str());//对应第一种方案定义的字母
				int y1 = atoi(str1.substr(8,7).c_str());//对应第二种方案定义的字母
				int i1 = BinarySearchChar(n1);//判断该氨基酸是否在20个标准氨基酸内
				if (i1==-1||x1==-999||y1==-999)
				{
					getline(datain,str1);
					continue;
				}//如果不是标准氨基酸，或者位于蛋白质分子序列的两端，读取文件下一行
				string str2;//记录后氨基酸的信息
				getline(datain,str2);
				if (str2.length()>10)//如果读取的不是pdb文件名
				{
					char n2 = str2.substr(0,1).c_str()[0];
					int x2 = atoi(str2.substr(3,5).c_str());
					int y2 = atoi(str2.substr(8,7).c_str());
					int i2 = BinarySearchChar(n2);
					if (i2==-1||x2==-999||y2==-999)
					{
						if(!datain.eof()) getline(datain,str1);
						continue;
					}//如果不是标准氨基酸，或者位于蛋白质分子序列的两端，读取文件下一行为前氨基酸
					//c1.add(i1,x1,i2,x2);
					c2.add(i1,y1,i2,y2);
					str1 = str2;//此时的后氨基酸变为下一次的前氨基酸
				}
			}
			//如果读取的是文件名，继续读取下一行
			else if(!datain.eof()) getline(datain,str1);
			cout<<"这是第："<<t1<<"行！"<<endl;
		}
		//c1.comput();
		c2.comput();
		//ofstream dataout1("G:\\task3\\result\\probability1.txt");
		//c1.print(dataout1);
		//ofstream dataout2("G:\\task3\\result\\probability2.txt");
		//c2.print(dataout2);
		initial();
		vector<HR_genel> allhr;
		ofstream of1("G:\\task3\\result\\probability12.txt");
		ofstream of2("G:\\task3\\result\\63.txt");
		prod_hrall(allhr,c2);
		rmsd_allhr(allhr);
		count6(allhr,of1,of2);
	}	
}

void main0(){
	counter c1(false);
	c1.comput();
}

