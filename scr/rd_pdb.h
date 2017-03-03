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
#define NOSEQNUM -9999	//���������ȡ�����ı��
#define NOCHAIN -1		//����ʶȡ�����ı��
//����һ���Ӵ�
typedef struct S_Contact
{
	int i1;
	int i2;
} Contact;
//������ά�ռ��е�һ����
class Point
{
public:
	Point();
public:
	double x,y,z;
};
//һ��ԭ��
class Atom
{
public:
	Atom();
	//������Ϣ
	Point pt;					//ԭ�ӵ�����
	string name;				//ԭ������
	int chain_seqm;				//ԭ�������е�˳���
	string residule_name;		//ԭ�����ڵİ���������
	char c_resname;				//������ĵ���ĸ����
	int residule_chain_seqm;	//ԭ�����ڵİ�����˳��ţ����ļ��ж���
	char chain_identifier;		//ԭ�����ڵ�����ʶ
	double occupancy;			//occupancy
	double tempFactor;			//tempFactor
	string segID;				//segID
	string element;				//element
	string charge;				//charge

	int residule_seqm;			//ԭ�����ڵİ��������ܰ��������е��±꣬��������и�ֵ����0��ʼ

	//������Ϣ
	double radium;				//ԭ�Ӱ뾶
};
//һ��������
class Residue
{
public:
	Residue();
	void Clear();
public:
	//������Ϣ
	string name;				//�����������
	char c_name;				//������ĵ���ĸ����	
	int chain_seqm;				//�����������е����
	int num_atom;				//������ԭ����
	vector<int> m_atoms;		//������ԭ�����
	int atom_index[NUM_ATOM_TYPE];	//��Ӧ��ԭ�����ܵ�ԭ���е���ţ���Դ�ļ�ʱ��ֵ��-1��ʾ�ð������в�������ԭ��

	//�������Ϣ
	double alpha,phi,psi,omega,chi_1,chi_2,chi_3,chi_4,chi_5;//alpha��ʾ�����ĸ�Caԭ�ӵ������

	//�Ӵ���Ϣ
	int num_contact;

	//������Ϣ
	double main_radium;			//�����뾶
	Point main_center;			//��������
	double side_radium;			//�����뾶
	Point side_center;			//��������
	

};
class Chain
{
public:
	Chain();
	void Clear();
	int AtomIndexInAA(int aa_num,string atom_name);//

	//contact���
	int ObtainContact(vector<Contact>& all_contact);	//������е�contact
	int SaveContact(const string& filename, const vector<Contact>& all_contact);	//��contact��Ϣд���ļ�

	//��������
	bool CalculateDihedral();//����ÿ������������е������
	bool WriteDihedral(const string filename);//���������Ϣд���ļ���
	bool Pre_Dihedral(vector<vector<double> >& preference);//ͳ�ư�����ԹǼ�Ťת�ǵĳ��ִ���

	//�Ӵ��������
	bool Cal_Num_Contact();//����ÿ��������ĽӴ�����
	bool Pre_Num_Contact(vector<vector<double> >& preference);//ͳ�ư�����ԽӴ��������͵ĳ��ִ���

	//ԭ���������
	bool Pre_Num_Atomtype1(vector<double >& preference);//�ۼ�ԭ������1�����ܵĸ�����ִ���
	bool Pre_Num_Atomtype2(vector<double >& preference);//�ۼ�ԭ������2�����ܵĸ�����ִ���

	//�����ؽ����
	bool Rebuildxyz_Dihedral();//��������Ƕ���������ؽ�
public:
	//������Ϣ
	char chain_identifier;		//����ʶ
	vector<Residue> m_residues;	//���еİ�����
	vector<Atom> m_atoms;		//���е�ԭ��
	int num_atom;				//ԭ�Ӹ���
	int num_residue;			//���������
	string seq;					//����������

	//������Ϣ

};
class Protein
{
public:
	Protein();
	void Clear();					//������е���
	bool LoadFromPDB(string filename);//��PDB�ļ��ж�����Ϣ
	bool SaveToFile(string filename);//��������Ϣ�����ļ���
	bool SaveChains(string dirname="AUTO",string filename="AUTO");//��ÿ��������Ϣ�ֱ�����ļ���
	bool WriteSeqFasta(string filename);//��������Ϣ��Fasta��ʽд���ļ���
	bool WriteSeqFasta_total(string filename);//��������������д��һ��

	int Get_Num_chain();//����������
	void Get_PDBID(string& id);//���PDBID
	bool GetSeq(int chain,string& seq);//��õ�i������������Ϣ

	//contact���
	int ObtainContact(int chainnumber,vector<Contact>& all_contact);	//��õ�i���������е�contact
	int SaveContact(int chainnumber,const string& filename, const vector<Contact >& all_contact);	//����i������contact��Ϣд���ļ�

	

private:
	
	void process_header(const string& line);
	void process_source(const string& line);
	void process_atom(const string& line,vector<Atom>& all_atom);
	void process_field(const string& line);//������ɢ���ֶΣ�������¼����
	bool AddAtom(const vector<Atom>& all_atom);
	void Atom2Line(const Atom& atom,string& line);
	void Atom2Ter(const Atom& atom,string& line);

	void OutputPreField(ofstream& file);//���atom�ֶ�ǰ���ֶ�
	void OutputSucField(ofstream& file);//���atom�ֶκ���ֶ�
	


public:
	//������Ϣ
	string PDB_id;
	vector<Chain> m_chains;		//���е���
	int num_chain;				//������

	vector<string> header;		//�洢���е�HEADER�ֶ�
	vector<string> title;		//�洢���е�TITLE�ֶ�
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


	//������Ϣ
};
bool IsDigit(string str);
#endif