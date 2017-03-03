#ifndef _RD_PDB_H
#define _RD_PDB_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <conio.h>
#include <math.h>
#include "amino_acid.h"
#include "direct.h"
#include "io.h"

using namespace std;
#define NOSEQNUM -9999	//氨基酸序号取不到的标记
#define NOCHAIN -1		//链标识取不到的标记
//代表一个接触
typedef struct S_Contact
{
	int i1;
	int i2;
} Contact;
//代表三维空间中的一个点
class Point
{
public:
	Point();
public:
	double x,y,z;
};
//一个原子
class Atom
{
public:
	Atom();
	//基本信息
	Point pt;					//原子的坐标
	string name;				//原子名称
	int chain_seqm;				//原子在链中的顺序号
	string residule_name;		//原子所在的氨基酸名称
	char c_resname;				//氨基酸的单字母名称
	int residule_chain_seqm;	//原子所在的氨基酸顺序号，从文件中读入
	char chain_identifier;		//原子所在的链标识
	double occupancy;			//occupancy
	double tempFactor;			//tempFactor
	string segID;				//segID
	string element;				//element
	string charge;				//charge

	int residule_seqm;			//原子所在的氨基酸在总氨基酸序列的下标，读入过程中赋值，从0开始

	//附加信息
	double radium;				//原子半径
};
//一个氨基酸
class Residue
{
public:
	Residue();
	void Clear();
public:
	//基本信息
	string name;				//氨基酸的名称
	char c_name;				//氨基酸的单字母名称	
	int chain_seqm;				//氨基酸在链中的序号
	int num_atom;				//包含的原子数
	vector<int> m_atoms;		//包含的原子序号
	int atom_index[NUM_ATOM_TYPE];	//相应的原子在总的原子中的序号，读源文件时赋值，-1表示该氨基酸中不包含该原子

	//两面角信息
	double alpha,phi,psi,omega,chi_1,chi_2,chi_3,chi_4,chi_5;//alpha表示连续四个Ca原子的两面角

	//接触信息
	int num_contact;

	//附加信息
	double main_radium;			//主链半径
	Point main_center;			//主链中心
	double side_radium;			//侧链半径
	Point side_center;			//侧链中心
	

};
class Chain
{
public:
	Chain();
	void Clear();
	int AtomIndexInAA(int aa_num,string atom_name);//

	//contact相关
	int ObtainContact(vector<Contact>& all_contact);	//获得所有的contact
	int SaveContact(const string& filename, const vector<Contact>& all_contact);	//将contact信息写入文件

	//两面角相关
	bool CalculateDihedral();//计算每个氨基酸的所有的两面角
	bool WriteDihedral(const string filename);//将两面角信息写入文件中
	bool Pre_Dihedral(vector<vector<double> >& preference);//统计氨基酸对骨架扭转角的出现次数

	//接触个数相关
	bool Cal_Num_Contact();//计算每个氨基酸的接触个数
	bool Pre_Num_Contact(vector<vector<double> >& preference);//统计氨基酸对接触个数类型的出现次数

	//原子势能相关
	bool Pre_Num_Atomtype1(vector<double >& preference);//累加原子类型1的势能的各项出现次数
	bool Pre_Num_Atomtype2(vector<double >& preference);//累加原子类型2的势能的各项出现次数

	//坐标重建相关
	bool Rebuildxyz_Dihedral();//根据两面角对坐标进行重建
public:
	//基本信息
	char chain_identifier;		//链标识
	vector<Residue> m_residues;	//所有的氨基酸
	vector<Atom> m_atoms;		//所有的原子
	int num_atom;				//原子个数
	int num_residue;			//氨基酸个数
	string seq;					//氨基酸序列

	//附加信息

};
class Protein
{
public:
	Protein();
	void Clear();					//清空所有的链
	bool LoadFromPDB(string filename);//从PDB文件中读入信息
	bool SaveToFile(string filename);//将坐标信息存在文件中
	bool SaveChains(string dirname="AUTO",string filename="AUTO");//将每个链的信息分别存入文件中
	bool WriteSeqFasta(string filename);//将序列信息以Fasta格式写入文件中
	bool WriteSeqFasta_total(string filename);//将所有链的序列写在一起

	int Get_Num_chain();//返回链总数
	void Get_PDBID(string& id);//获得PDBID
	bool GetSeq(int chain,string& seq);//获得第i个链的序列信息

	//contact相关
	int ObtainContact(int chainnumber,vector<Contact>& all_contact);	//获得第i个链的所有的contact
	int SaveContact(int chainnumber,const string& filename, const vector<Contact >& all_contact);	//将第i个链的contact信息写入文件

	

private:
	
	void process_header(const string& line);
	void process_source(const string& line);
	void process_atom(const string& line,vector<Atom>& all_atom);
	void process_field(const string& line);//处理零散的字段，仅仅记录数据
	bool AddAtom(const vector<Atom>& all_atom);
	void Atom2Line(const Atom& atom,string& line);
	void Atom2Ter(const Atom& atom,string& line);

	void OutputPreField(ofstream& file);//输出atom字段前的字段
	void OutputSucField(ofstream& file);//输出atom字段后的字段
	


public:
	//基本信息
	string PDB_id;
	vector<Chain> m_chains;		//所有的链
	int num_chain;				//链总数

	vector<string> header;		//存储所有的HEADER字段
	vector<string> title;		//存储所有的TITLE字段
	vector<string> compnd;
	vector<string> source;
	vector<string> author;
	vector<string> revdat;
	vector<string> jrnl;
	vector<string> remark;
	vector<string> dbref;
	vector<string> seqadv;
	vector<string> seqres;
	vector<string> het;
	vector<string> hetnam;
	vector<string> formul;
	vector<string> helix;
	vector<string> sheet;
	vector<string> link;
	vector<string> cispep;
	vector<string> origx1;
	vector<string> origx2;
	vector<string> origx3;
	vector<string> scale1;
	vector<string> scale2;
	vector<string> scale3;
	vector<string> hetatm;
	vector<string> concet;
	vector<string> master;
	vector<string> anisou;
	vector<string> cryst1;
	vector<string> expdta;
	vector<string> ftnote;
	vector<string> hetsyn;
	vector<string> hydbnd;
	vector<string> jnrl;
	vector<string> keywds;
	vector<string> modres;
	vector<string> mtrix1;
	vector<string> mtrix2;
	vector<string> mtrix3;	
	vector<string> site;
	vector<string> sltbrg;
	vector<string> sprsde;
	vector<string> ssbond;
	vector<string> turn;
	vector<string> tvect;
	vector<string> conect;
	vector<string> end;


	//附加信息
};
bool IsDigit(string str);
#endif