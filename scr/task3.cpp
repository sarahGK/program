// task3.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "globfuc.h"
#include "container.h"
#include "hr_Info.h"
#include "readin_Array.h"


void  main()
{
	ifstream datain;
	datain.open("E:\\talk\\biology\\task2\\alphabet.txt",ios::in);//��hr����ĸ���ļ�
	if(!datain.is_open()){
		cout<<"can't open the file: E:\\talk\\biology\\task2\\alphabet.txt"<<endl;
	}else{
		int t1 = 0;
		string str1;	//��¼ǰ���������Ϣ
		counter c2(false);//����ڶ�����ĸ���巽���ļ�����
		cout<<"��ʼ��ȡ�ļ�������"<<endl;
		getline(datain,str1);//��ȡ�ļ��ĵ�һ��
		while(!datain.eof()){
			t1 ++;
			if(str1.length()>10){//�����ȡ�Ĳ���pdb�ļ���
				char n1 = str1.substr(0,1).c_str()[0];//���������ĸ��ʾ
				int x1 = atoi(str1.substr(3,5).c_str());//��Ӧ��һ�ַ����������ĸ
				int y1 = atoi(str1.substr(8,7).c_str());//��Ӧ�ڶ��ַ����������ĸ
				int i1 = BinarySearchChar(n1);//�жϸð������Ƿ���20����׼��������
				if (i1==-1||x1==-999||y1==-999)
				{
					getline(datain,str1);
					continue;
				}//������Ǳ�׼�����ᣬ����λ�ڵ����ʷ������е����ˣ���ȡ�ļ���һ��
				string str2;//��¼�󰱻������Ϣ
				getline(datain,str2);
				if (str2.length()>10)//�����ȡ�Ĳ���pdb�ļ���
				{
					char n2 = str2.substr(0,1).c_str()[0];
					int x2 = atoi(str2.substr(3,5).c_str());
					int y2 = atoi(str2.substr(8,7).c_str());
					int i2 = BinarySearchChar(n2);
					if (i2==-1||x2==-999||y2==-999)
					{
						if(!datain.eof()) getline(datain,str1);
						continue;
					}//������Ǳ�׼�����ᣬ����λ�ڵ����ʷ������е����ˣ���ȡ�ļ���һ��Ϊǰ������
					//c1.add(i1,x1,i2,x2);
					c2.add(i1,y1,i2,y2);
					str1 = str2;//��ʱ�ĺ󰱻����Ϊ��һ�ε�ǰ������
				}
			}
			//�����ȡ�����ļ�����������ȡ��һ��
			else if(!datain.eof()) getline(datain,str1);
			cout<<"���ǵڣ�"<<t1<<"�У�"<<endl;
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

