#include "stdafx.h"
#include "rdFasta.h"

/*******************************************************
 * ��Fasta�ļ��ж�������
 * ���ɹ��Ķ�����һ�����У�����true
 * ���򷵻�false,��ʾ���ļ��в���������
 ******************************************************/
bool ReadFastaSeq(ifstream& file,string &seq)
{
	string line;
	char c;
	line.push_back('\0');
	if(!file.is_open())
		return false;
	seq.clear();
	while(((int)line.size()>0)&&(line[0]!='>')&&(!file.eof()))
	{
		getline(file,line);
	}
	if(file.eof())
		return false;
	file.get(c);
	while((c!='>')&&(!file.eof()))
	{
		file.unget();
		getline(file,line);seq+=line;		
		file.get(c);
	}
	if(!file.eof())
		file.unget();
	return true;	
}
//��������Ϣһ�������
bool ReadFastaSeq(ifstream& file,string &seq,string& miaoshu)
{
	string line;
	char c;
	line.push_back('\0');
	if(!file.is_open())
		return false;
	seq.clear();	
	while(((int)line.size()>0)&&(line[0]!='>')&&(!file.eof()))
	{
		getline(file,line);
	}
	if(file.eof())
		return false;
	if(line.size()>0&&line[0]=='>')
	{
		line.erase(0,1);
		miaoshu=line;
	}
	file.get(c);
	while((c!='>')&&(!file.eof()))
	{
		file.unget();
		getline(file,line);seq+=line;		
		file.get(c);
	}
	if(!file.eof())
		file.unget();
	return true;	
}
/***********************************************
* ������ʾһ�������, ����ʶ��Ϊ������
* ��Ϣ�ĵ�5����ĸ,����>1ACXB
***********************************************/
bool ReadFastaSeq(ifstream& file,string &seq,char& cc)
{
	string line;
	char c;
	line.push_back('\0');
	
	if(!file.is_open())
		return false;
	seq.clear();	
	while(((int)line.size()>0)&&(line[0]!='>')&&(!file.eof()))
	{
		getline(file,line);
	}
	if(file.eof())
		return false;
	if((int)line.length()>=6)
		cc=line[5];
	file.get(c);
	while((c!='>')&&(!file.eof()))
	{
		file.unget();
		getline(file,line);seq+=line;		
		file.get(c);
	}
	if(!file.eof())
		file.unget();
	return true;	
}
/***********************************************
* ��������Ϣ������ʾһ�������
***********************************************/
bool ReadFastaSeq(ifstream& file,string &seq,string& miaoshu,char& cc)
{
	string line;
	char c;
	line.push_back('\0');
	
	if(!file.is_open())
		return false;
	seq.clear();	
	miaoshu.clear();
	while(((int)line.size()>0)&&(line[0]!='>')&&(!file.eof()))
	{
		getline(file,line);
	}
	if(file.eof())
		return false;	
	if((int)line.length()>=6)
		cc=line[5];
	if(line.size()>0&&line[0]=='>')
	{
		line.erase(0,1);
		miaoshu=line;
	}
	file.get(c);
	while((c!='>')&&(!file.eof()))
	{
		file.unget();
		getline(file,line);seq+=line;		
		file.get(c);
	}
	if(!file.eof())
		file.unget();
	return true;	
}
/******************************************
 * ���ļ���д��ָ��������
 * �Ḳ��ԭ�����ļ�
 * ������ϢΪ�ļ���
 * �����ʽΪÿ��70��������
 *******************************************/
bool WriteFastaSeq(string& filename,string &seq)
{
	ofstream file;
	int i,j,end,len;	
	if(filename.empty())
		return false;
	file.open(filename.c_str());
	if(!file.is_open())
		return false;
	file<<'>'<<filename<<endl;
	int row,col;
	len=(int)seq.length();
	row=len/70+1;
	col=len%70;
	for(i=0;i<row;i++)
	{
		if(i==(row-1))
			end=col;
		else
			end=70;
		for(j=0;j<end;j++)
			file<<seq[i*70+j];
		file<<endl;		
	}
	return true;
}
bool WriteFastaSeq(string& filename,string &seq,string& miaoshu)
{
	ofstream file;
	int i,j,end,len;	
	if(filename.empty())
		return false;
	file.open(filename.c_str());
	if(!file.is_open())
		return false;
	if(miaoshu.empty())
		file<<'>'<<filename<<endl;
	else
	{
		if(miaoshu[0]=='>')
			file<<miaoshu<<endl;
		else
			file<<">"<<miaoshu<<endl;
	}
	int row,col;
	len=(int)seq.length();
	row=len/70+1;
	col=len%70;
	for(i=0;i<row;i++)
	{
		if(i==(row-1))
			end=col;
		else
			end=70;
		for(j=0;j<end;j++)
			file<<seq[i*70+j];
		file<<endl;		
	}
	return true;
}
/**************************************
 * ���ļ���׷��ָ��������
 * �ļ����봦�ڴ�״̬
 *************************************/
bool WriteFastaSeq(ofstream& file,string &seq,string& miaoshu)
{
	int i,j,end,len;
	
	if(!file.is_open())
		return false;
	int row,col;
	if(miaoshu.empty())
		file<<'>'<<endl;
	else
	{
		if(miaoshu[0]=='>')
			file<<miaoshu<<endl;
		else
			file<<">"<<miaoshu<<endl;
	}
	len=(int)seq.length();
	row=len/70+1;
	col=len%70;
	for(i=0;i<row;i++)
	{
		if(i==(row-1))
			end=col;
		else
			end=70;
		for(j=0;j<end;j++)
			file<<seq[i*70+j];
		file<<endl;		
	}
	return true;
}
bool WriteFastaSeq(ofstream& file,string &seq)
{
	int i,j,end,len;
	
	
	if(!file.is_open())
		return false;
	int row,col;
	
	len=(int)seq.length();
	row=len/70+1;
	col=len%70;
	for(i=0;i<row;i++)
	{
		if(i==(row-1))
			end=col;
		else
			end=70;
		for(j=0;j<end;j++)
			file<<seq[i*70+j];
		file<<endl;		
	}
	return true;
}




	

