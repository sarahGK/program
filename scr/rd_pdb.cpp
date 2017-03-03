#include "stdafx.h"
#include "rd_pdb.h"
#include "amino_acid.h"
#include "vector.h"
#include "rdfasta.h"
Point::Point()
{
	x=0;y=0;z=0;
}
Atom::Atom()
{
	chain_seqm=0;
	c_resname='x';
	residule_chain_seqm=0;
	chain_identifier=-1;
	radium=0;
}
Residue::Residue()
{
	c_name='X';
	chain_seqm=0;
	num_atom=0;	
	alpha=phi=psi=omega=chi_1=chi_2=chi_3=chi_4=chi_5=-999.0;
	for(int i=0;i<NUM_ATOM_TYPE;i++)
		atom_index[i]=-1;
	num_contact=0;
}
Chain::Chain()
{
	chain_identifier='X';
	num_residue=0;
	num_atom=0;
}
Protein::Protein()
{
	num_chain=0;
}
void Residue::Clear()
{
	name.clear();
	c_name='X';
	chain_seqm=0;
	num_atom=0;
	m_atoms.clear();
	alpha=phi=psi=omega=chi_1=chi_2=chi_3=chi_4=chi_5=-999.0;
	for(int i=0;i<NUM_ATOM_TYPE;i++)
		atom_index[i]=-1;
	num_contact=0;
}
void Chain::Clear()
{
	chain_identifier=-1;
	m_residues.clear();
	num_residue=0;
	m_atoms.clear();
	num_atom=0;
	seq.clear();
	
}
void Protein::Clear()
{
	m_chains.clear();
	num_chain=0;
	PDB_id.clear();
}
/*************************************************
 * ��PDB�ļ��ж�����Ϣ��
 * �ɹ�ʱ����true����ʱ����false 
 *************************************************/
bool  Protein::LoadFromPDB(string filename)
{
	ifstream filein;
	string line;
	bool model=false;//�Ƿ��Ƕ��ģ�ͣ��������һ��

	vector<Atom> all_atom;
	if(filename.empty())
		return false;
	filein.open(filename.c_str(),ios::in);
	if(!filein.is_open())
		return false;

	Clear();
	while(!filein.eof())
	{
		getline(filein,line);		
		/*if(strcmp(line.substr(0,6).c_str(),"HEADER")==0) 
			process_header(line);
		else if(strcmp(line.substr(0,6).c_str(),"SOURCE")==0) 
			process_source(line);*/
		 if(strcmp(line.substr(0,6).c_str(),"ATOM  ")==0) 
		{
			if(!model)
				process_atom(line,all_atom);
		}
		else if(strcmp(line.substr(0,6).c_str(),"ENDMDL")==0)
			model=true;
		else process_field(line);
	}

	AddAtom(all_atom);
	return true;
}
/*********************************************************
 * ����Ϣ�����ļ���,�ɹ�ʱ����true
 ********************************************************/
bool Protein::SaveToFile(string filename)
{
	ofstream file;
	string line;
	Atom atom;
	int i,j;
	if(filename.empty())
		return false;
	file.open(filename.c_str(),ios::out);
	if(!file.is_open())
		return false;
	/*line="HEADER";
	line.resize(80);
	for(i=0;i<(int)PDB_id.length();i++)
		line[62+i]=PDB_id[i];
	file<<line<<endl;*/
	
	OutputPreField(file);
	int num_atom;
	for(i=0;i<num_chain;i++)
	{
		num_atom=m_chains[i].num_atom;
		for(j=0;j<num_atom;j++)
		{
			atom=m_chains[i].m_atoms[j];
			Atom2Line(atom,line);
			file<<line<<endl;
		}
		Atom2Ter(atom,line);
		file<<line<<endl;
	}
	OutputSucField(file);
	return true;
}
	



void Protein::process_header(const string& line)
{
	PDB_id=line.substr(62,4);
	process_field(line);
}
void Protein::process_source (const string& line)
{
	process_field(line);
}
void Protein::process_atom (const string& line, vector<Atom>& all_atom)
{
	int i,len;
	len=(int)line.size();
	Atom atom;
	string subline;
	if(len<7)
		return;
	subline=line.substr(6,5);
	atom.chain_seqm=atoi(subline.c_str());
	subline=line.substr(12,4);
	for(i=0;i<4;i++)
	{
		if(subline[i]!=' ')
			atom.name.push_back(subline[i]);
		Upper(atom.name);
	}
	subline=line.substr(17,3);
	atom.residule_name=subline;
	Upper(atom.residule_name);
	atom.c_resname=residulename321(subline.c_str());
	atom.chain_identifier=line.length()>=22?line[21]:' ';
	subline=line.substr(22,5);
	if(IsDigit(subline))
		atom.residule_chain_seqm=atoi(subline.c_str());
	else
		atom.residule_chain_seqm=-1;
	subline=line.substr(30,8);
	atom.pt.x=atof(subline.c_str());
	subline=line.substr(38,8);
	atom.pt.y=atof(subline.c_str());
	subline=line.substr(46,8);
	atom.pt.z=atof(subline.c_str());
	subline=line.substr(54,6);
	atom.occupancy=atof(subline.c_str());
	subline=line.substr(60,6);
	atom.tempFactor=atof(subline.c_str());
	if(len>=73)
		atom.segID=line.substr(72,4);
	if(len>=77)
		atom.element=line.substr(76,2);
	if(len>=79)
		atom.charge=line.substr(78,2);
	all_atom.push_back(atom);
}
/************************************************************
 * �Ѷ����ԭ����ӵ��������У�ע��ԭ�ӱ����ǰ���˳������
 ************************************************************/
bool Protein::AddAtom(const vector<Atom>& all_atom)
{
	int len,i,index;
	Atom atom;
	Chain chain;
	Residue residue;

	chain.chain_identifier=NOCHAIN;
	residue.chain_seqm=NOSEQNUM;
	len=(int)all_atom.size();
	if(len<=0)//���û��Ժ�ӣ�ֱ�ӷ���
		return true;
	for(i=0;i<len;i++)
	{
		atom=all_atom[i];
		if(atom.residule_chain_seqm!=residue.chain_seqm)
		{
			if(residue.chain_seqm!=NOSEQNUM)
			{
				chain.m_residues.push_back(residue);
				chain.num_residue++;			
				chain.seq.push_back(residue.c_name);
				residue.Clear();
			}
				residue.chain_seqm=atom.residule_chain_seqm;
				residue.name=atom.residule_name;
				residue.c_name=atom.c_resname;
				residue.num_atom=0;			
		}
		if(atom.chain_identifier!=chain.chain_identifier)
		{
			if(chain.chain_identifier!=NOCHAIN)
			{
				m_chains.push_back(chain);
				chain.Clear();			
				num_chain++;
			}
				chain.chain_identifier=atom.chain_identifier;
				chain.num_residue=0;	
				chain.num_atom=0;			
		}		
		atom.residule_seqm=chain.num_residue;
		chain.m_atoms.push_back(atom);
		residue.m_atoms.push_back(chain.num_atom);
		//residue.chain_seqm=atom.residule_chain_seqm;
		index=Atomname2ID(atom.name.c_str());
		if(index!=-1)
		{
			residue.atom_index[index]=chain.num_atom;
		}
		chain.num_atom++;		
		residue.num_atom++;
	}
	chain.m_residues.push_back(residue);
	chain.num_residue++;
	chain.seq.push_back(residue.c_name);
	m_chains.push_back(chain);
	num_chain++;
	return true;
}

void Protein::Atom2Line(const Atom& atom,string& line)
{
	char str[81];
	char name[5]="    ";
	if((int)atom.name.length()==4)//���ԭ��������Ϊ4
	{
		for(int i=0;i<(int)atom.name.length();i++)
			name[i]=atom.name[i];
	}
	else//������Ȳ�Ϊ4�����һ��Ϊ�ո�
	{
		for(int i=0;i<(int)atom.name.length()&&i<3;i++)
			name[i+1]=atom.name[i];
	}
	memset(str,' ',81);
	str[80]='\0';
	sprintf(str,"%s%5d %s %s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %s%s%s","ATOM  ",atom.chain_seqm,name,atom.residule_name.c_str(),atom.chain_identifier,atom.residule_chain_seqm,atom.pt.x,atom.pt.y,atom.pt.z,atom.occupancy,atom.tempFactor,atom.segID.c_str(),atom.element.c_str(),atom.charge.c_str());
	line=str;
}
void Protein::Atom2Ter(const Atom& atom,string& line)
{
	char str[81];	
	memset(str,' ',81);
	str[80]='\0';
	sprintf(str,"%s%5d      %s %c%4d                                                      ","TER   ",atom.chain_seqm+1,atom.residule_name.c_str(),atom.chain_identifier,atom.residule_chain_seqm);
	line=str;
}
/********************************************
  * ��õ�i���������е�contact,����contact����
 ********************************************/
int Protein::ObtainContact(int chainnumber,vector<Contact>& all_contact)
{
	if(chainnumber<num_chain)
	{
		return m_chains[chainnumber].ObtainContact(all_contact);
	}
	else
		return 0;
}
/**************************************************
 * ����i������contact��Ϣ�����ļ���
 * �ɹ�ʱ����0
 * ��������ʱ������-1
 * ***********************************************/
int Protein::SaveContact(int chainnumber,const string& filename,const vector<Contact>& all_contact)
{
	if(chainnumber<num_chain)
		return m_chains[chainnumber].SaveContact(filename,all_contact);
	else
		return -1;
}
int Chain::ObtainContact(vector<Contact >& all_contact)
{
	int i,j,len;
	int id;
	double sider,d,t;
	bool findca,findcb;
	Point capt,cbpt,sidept;
	len=num_residue;
	//��һ��ɨ�裬���ÿ��������Ĳ����뾶�����ģ�����������ca��cb�����������뾶�������ڸʰ����Ƚ���������Ca��
	for(i=0;i<len;i++)
	{
		id=Rname2ID(m_residues[i].c_name);
		sider=AA_SideR[id];
		m_residues[i].side_radium=sider;
		findca=findcb=false;
		for(j=0;j<m_residues[i].num_atom;j++)
		{
			if(m_atoms[m_residues[i].m_atoms[j]].name=="CA")
			{
				capt=m_atoms[m_residues[i].m_atoms[j]].pt;
				findca=true;
			}
			if(m_atoms[m_residues[i].m_atoms[j]].name=="CB")
			{
				cbpt=m_atoms[m_residues[i].m_atoms[j]].pt;
				findcb=true;
			}
		}
		if(findca&&findcb)
		{
			d=Distance(capt,cbpt);
			t=sider/d;
			sidept.x=capt.x+(cbpt.x-capt.x)*t;
			sidept.y=capt.y+(cbpt.y-capt.y)*t;
			sidept.z=capt.z+(cbpt.z-capt.z)*t;
			m_residues[i].side_center=sidept;
		}
		else if(findca)
		{
			m_residues[i].side_center=capt;
		}				
	}
	Contact contact;
	double r1,r2;
	for(i=0;i<len;i++)
	{
		r1=m_residues[i].side_radium;
		for(j=i+1;j<len;j++)
		{
			r2=m_residues[j].side_radium;
			d=Distance(m_residues[i].side_center,m_residues[j].side_center);
			if(d<r1+r2+2.8)
			{
				contact.i1=i;
				contact.i2=j;
				all_contact.push_back(contact);
			}
		}
	}
	return (int)all_contact.size();
}
/*******************************************
 * ��һ������contact��Ϣ�����ļ���
 * �洢��ʽΪ��
 * ����������
 * contact����
 * ���е�contact
 * �ɹ�ʱ����0�����򣬷��ط�0����ʾ
 * 1 �򲻿��ļ�
 * 2 �޷�д���ļ�
 ******************************************/
int Chain::SaveContact(const string& filename,const vector<Contact>& all_contact)
{
	ofstream fileout;
	int len;
	if(filename.empty())
		return 1;
	fileout.open(filename.c_str(),ios::out);
	if(!fileout.is_open())
		return 1;
	fileout<<seq.c_str()<<endl;
	len=(int)all_contact.size();
	fileout<<len<<endl;
	for(int i=0;i<len;i++)
		fileout<<all_contact[i].i1<<' '<<all_contact[i].i2<<' ';	
	return 0;
}
int Protein::Get_Num_chain()
{
	return num_chain;
}
void Protein::Get_PDBID(string& id)
{
	id=PDB_id;
}
bool Protein::GetSeq(int chain,string& seq)
{
	if(chain>num_chain)
		return false;
	seq=m_chains[chain].seq;
	return true;
}
bool Chain::CalculateDihedral()
{
	int len,i,index,atom_index;
	Point p1,p2,p3,p4,p5,p6;
	len=(int)m_residues.size();
	for(i=0;i<len;i++)
	{
		//����alpha
		if(i==1)
		{
			index=Atomname2ID("CA");
			atom_index=m_residues[i-1].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"���������ʱ�����ڰ�����"<<i-1<<"���Ҳ���CAԭ��"<<endl;
				return false;
			}
			p1=m_atoms[atom_index].pt;
			atom_index=m_residues[i].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"���������ʱ�����ڰ�����"<<i<<"���Ҳ���CAԭ��"<<endl;
				return false;
			}
			p2=m_atoms[atom_index].pt;
			atom_index=m_residues[i+1].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"���������ʱ�����ڰ�����"<<i+1<<"���Ҳ���CAԭ��"<<endl;
				return false;
			}
			p3=m_atoms[atom_index].pt;
			atom_index=m_residues[i+2].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"���������ʱ�����ڰ�����"<<i+2<<"���Ҳ���CAԭ��"<<endl;
				return false;
			}
			p4=m_atoms[atom_index].pt;
			m_residues[i].alpha=Dihedralangle(p1,p2,p3,p4);
		}
		else if((i>1)&&(i<len-2))
		{
			p1=p2;
			p2=p3;
			p3=p4;
			index=Atomname2ID("CA");
			atom_index=m_residues[i+2].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"���������ʱ�����ڰ�����"<<i+2<<"���Ҳ���CAԭ��"<<endl;
				return false;
			}
			p4=m_atoms[atom_index].pt;
			m_residues[i].alpha=Dihedralangle(p1,p2,p3,p4);
		}

		//����phi,psi,omega
		//if(i>0)
		{
			if(i>0)
			{
				index=Atomname2ID("C");
				atom_index=m_residues[i-1].atom_index[index];
				if(atom_index==-1)
				{
					cout<<"���������ʱ�����ڰ�����"<<i-1<<"���Ҳ���Cԭ��"<<endl;
					return false;
				}
				p1=m_atoms[atom_index].pt;
			}
			index=Atomname2ID("N");
			atom_index=m_residues[i].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"���������ʱ�����ڰ�����"<<i<<"���Ҳ���Nԭ��"<<endl;
				return false;
			}
			p2=m_atoms[atom_index].pt;
			index=Atomname2ID("CA");
			atom_index=m_residues[i].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"���������ʱ�����ڰ�����"<<i<<"���Ҳ���CAԭ��"<<endl;
				return false;
			}
			p3=m_atoms[atom_index].pt;
			index=Atomname2ID("C");
			atom_index=m_residues[i].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"���������ʱ�����ڰ�����"<<i<<"���Ҳ���Cԭ��"<<endl;
				return false;
			}
			p4=m_atoms[atom_index].pt;
			if(i>0)
				m_residues[i].phi =Dihedralangle(p1,p2,p3,p4);
			if(i<len-1)
			{
				index=Atomname2ID("N");
				atom_index=m_residues[i+1].atom_index[index];
				if(atom_index==-1)
				{
					cout<<"���������ʱ�����ڰ�����"<<i+1<<"���Ҳ���Cԭ��"<<endl;
					return false;
				}
				p5=m_atoms[atom_index].pt;
				m_residues[i].psi=Dihedralangle(p2,p3,p4,p5);
			
				index=Atomname2ID("CA");
				atom_index=m_residues[i+1].atom_index[index];
				if(atom_index==-1)
				{
					cout<<"���������ʱ�����ڰ�����"<<i+1<<"���Ҳ���CAԭ��"<<endl;
					return false;
				}
				p6=m_atoms[atom_index].pt;
				m_residues[i].omega=Dihedralangle(p3,p4,p5,p6);
			}
		}
		//����X1-X5
	
	}
		return true;
}
bool Chain::WriteDihedral(const string filename)
{
	ofstream file(filename.c_str());
	if(!file.is_open())
		return false;
	int i,len;
	len=(int)m_residues.size();
	file<<"�����"<<endl;
	
	file<<"���"<<'\t'<<"����������"<<'\t'<<"alpha"<<'\t'<<"phi"<<'\t'<<"psi\t"<<"omega\t"<<endl;
	for(i=0;i<len;i++)
	{
		file<<i<<'\t'<<m_residues[i].name<<'\t';
		file.precision(7);
		file<<m_residues[i].alpha<<'\t';
		file.precision(7);
		file<<m_residues[i].phi<<'\t';
		file.precision(7);
		file<<m_residues[i].psi<<'\t';
		file.precision(7);
		file<<m_residues[i].omega <<'\t';

		file<<endl;
	}
	return true;
}
/*********************************************************
 * ��������ǽ��е�������ά�ṹ�������ؽ�
 * ���巽����ο�������Ϣѧѧϰ�ʼ�2�е�ר��8.8
 * ��Ҫ˵������,�������ؽ�������ʹ�õļ�������
 * �������������ʱ�ľ���ֵ���ؽ�������겢����
 * ��ԭ��������ȫ�غ�
 * У׼�ļ�������ȡ��AMBER���ܺ������������ʽṹ
 * ׼��һ�ġ�
 * ���⣬��ʹ�ǰ���PDB�ļ��и������������������
 * �ͼ��ǣ�Ҳ���ܺͱ�׼״̬�µ�ֵ����ȫһ��
 * Ŀǰ����ɶ�����ԭ��(����Oԭ��)�����ؽ�
 * ʵ������ʾ�������ؽ�������ʵ��ֵ�����Ǻ�
 * ���ǻ�����ۼ��������ԭ�ӵ�����ƫ��ϴ�
 * ��ԼΪ2��
 *********************************************************/
bool Chain::Rebuildxyz_Dihedral()
{
	int len=(int)m_residues.size();
	int i;
	int index,atom_index;
	Point qN,qCA,qC;//ǰ�氱�����Ca,C,Nԭ������
	Point N,CA,C,O;//��ǰ�������N,Ca,C,Oԭ������
	double theta,tao,bond_len;

	//�Ե�һ�����������Ӧԭ�Ӹ�ֵ��������ԭ������

	for(i=1;i<len-1;i++)//���ζ����İ������е�ԭ�ӽ��������ؽ�
	{
		//���㵱ǰ�������N
		index=Atomname2ID("N");
		atom_index=m_residues[i-1].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"�����ؽ�ʱʱ�����ڰ�����"<<i-1<<"���Ҳ���Nԭ��"<<endl;
			return false;
		}
		qN=m_atoms[atom_index].pt;

		index=Atomname2ID("CA");
		atom_index=m_residues[i-1].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"�����ؽ�ʱʱ�����ڰ�����"<<i-1<<"���Ҳ���CAԭ��"<<endl;
			return false;
		}
		qCA=m_atoms[atom_index].pt;
		index=Atomname2ID("C");
		atom_index=m_residues[i-1].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"�����ؽ�ʱʱ�����ڰ�����"<<i-1<<"���Ҳ���Cԭ��"<<endl;
			return false;
		}
		qC=m_atoms[atom_index].pt;
		
		//���㵱ǰ��Nԭ������
		theta=m_residues[i-1].psi;//�����
		tao=116.6;//Ca-C-N����
		bond_len=1.335;//C-N����
		N=LastXYZ(qN,qCA,qC,theta,tao,bond_len);
		index=Atomname2ID("N");
		atom_index=m_residues[i].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"�����ؽ�ʱʱ�����ڰ�����"<<i<<"���Ҳ���Nԭ��"<<endl;
			return false;
		}
		m_atoms[atom_index].pt=N;

		//�Ե�ǰ��CAԭ�ӽ��������ؽ�
		theta=m_residues[i-1].omega;
		tao=121.9;//C-N-CA����
		bond_len=1.449;//N-CA����
		CA=LastXYZ(qCA,qC,N,theta,tao,bond_len);
		index=Atomname2ID("CA");
		atom_index=m_residues[i].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"�����ؽ�ʱʱ�����ڰ�����"<<i<<"���Ҳ���CAԭ��"<<endl;
			return false;
		}
		m_atoms[atom_index].pt=CA;

		//�Ե�ǰ��Cԭ�ӽ��������ؽ�
		theta=m_residues[i].phi;
		tao=110.3;//N-CA-C����
		bond_len=1.522;//CA-C����
		C=LastXYZ(qC,N,CA,theta,tao,bond_len);
		index=Atomname2ID("C");
		atom_index=m_residues[i].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"�����ؽ�ʱʱ�����ڰ�����"<<i<<"���Ҳ���Cԭ��"<<endl;
			return false;
		}
		m_atoms[atom_index].pt=C;


	}
	return true;
}
bool Protein::WriteSeqFasta(string filename)
{
	if(filename.empty())
		return false;
	ofstream file;
	file.open(filename.c_str());
	if(!file.is_open())
		return false;
	string header,header1;
	header1=this->PDB_id;
	if(header1.empty())
		header1=filename;
	int i;
	for(i=0;i<(int)m_chains.size();i++)
	{
		header=header1+m_chains[i].chain_identifier;
		if(!WriteFastaSeq(file,m_chains[i].seq,header))
			return false;
	}
	return true;
}
bool Protein::WriteSeqFasta_total(string filename)
{
	if(filename.empty())
		return false;
	string seq;
	ofstream file;
	file.open(filename.c_str());
	if(!file.is_open())
		return false;
	string header,header1;
	header1=this->PDB_id;
	if(header1.empty())
		header1=filename;
	int i;
	for(i=0;i<(int)m_chains.size();i++)
	{		
		seq+=m_chains[i].seq;			
	}
	if(!WriteFastaSeq(file,seq,header1))
		return false;
	return true;
}
//������ɢ���ֶΣ�������¼����
void Protein::process_field(const string& line)
{
	string subline=line.substr(0,6);
	if(subline=="ANISOU")
		anisou.push_back(line);
	else if(subline=="AUTHOR")
		author.push_back(line);
	else if(subline=="CISPEP")
		cispep.push_back(line);
	else if(subline=="COMPND")
		compnd.push_back(line);
	else if(subline=="CONECT")
		conect.push_back(line);
	else if(subline=="CRYST1")
		cryst1.push_back(line);
	else if(subline=="DBREF ")
		dbref.push_back(line);
	else if(subline=="EXPDTA")
		expdta.push_back(line);
	else if(subline=="FORMUL")
		formul.push_back(line);
	else if(subline=="FTNOTE")
		ftnote.push_back(line);
	else if(subline=="HEADER")
		header.push_back(line);
	else if(subline=="HELIX ")
		helix.push_back(line);
	else if(subline=="HET   ")
		het.push_back(line);
	else if(subline=="HETATM")
		hetatm.push_back(line);
	else if(subline=="HETNAM")
		hetnam.push_back(line);
	else if(subline=="HETSYN")
		hetsyn.push_back(line);
	else if(subline=="HYDBND")
		hydbnd.push_back(line);
	else if(subline=="JNRL  ")
		jnrl.push_back(line);
	else if(subline=="JRNL  ")
		jrnl.push_back(line);
	else if(subline=="KEYWDS")
		keywds.push_back(line);
	else if(subline=="LINK  ")
		link.push_back(line);
	else if(subline=="MASTER")
		master.push_back(line);
	else if(subline=="MODRES")
		modres.push_back(line);
	else if(subline=="MTRIX1")
		mtrix1.push_back(line);
	else if(subline=="MTRIX2")
		mtrix2.push_back(line);
	else if(subline=="MTRIX3")
		mtrix3.push_back(line);
	else if(subline=="ORIGX1")
		origx1.push_back(line);
	else if(subline=="ORIGX2")
		origx2.push_back(line);
	else if(subline=="ORIGX3")
		origx3.push_back(line);
	else if(subline=="REMARK")
		remark.push_back(line);
	else if(subline=="REVDAT")
		revdat.push_back(line);
	else if(subline=="SCALE1")
		scale1.push_back(line);
	else if(subline=="SCALE2")
		scale2.push_back(line);
	else if(subline=="SCALE3")
		scale3.push_back(line);
	else if(subline=="SEQADV")
		seqadv.push_back(line);
	else if(subline=="SEQRES")
		seqres.push_back(line);
	else if(subline=="SHEET ")
		sheet.push_back(line);
	else if(subline=="SITE  ")
		site.push_back(line);
	else if(subline=="SLTBRG")
		sltbrg.push_back(line);
	else if(subline=="SOURCE")
		source.push_back(line);
	else if(subline=="SPRSDE")
		sprsde.push_back(line);
	else if(subline=="SSBOND")
		ssbond.push_back(line);
	else if(subline=="TITLE ")
		title.push_back(line);
	else if(subline=="TURN  ")
		turn.push_back(line);
	else if(subline=="TVECT")
		tvect.push_back(line);
	else if(subline=="END   ")
		end.push_back(line);
}
void Protein::OutputPreField(ofstream& file)
{
	int i;
	if(!header.empty())
	{
		for(i=0;i<(int)header.size();i++)
			file<<header[i]<<endl;
	}
	if(!title.empty())
	{
		for(i=0;i<(int)title.size();i++)
			file<<title[i]<<endl;
	}
	if(!compnd.empty())
	{
		for(i=0;i<(int)compnd.size();i++)
			file<<compnd[i]<<endl;
	}
	if(!source.empty())
	{
		for(i=0;i<(int)source.size();i++)
			file<<source[i]<<endl;
	}
	if(!keywds.empty())
	{
		for(i=0;i<(int)keywds.size();i++)
			file<<keywds[i]<<endl;
	}
	if(!expdta.empty())
	{
		for(i=0;i<(int)expdta.size();i++)
			file<<expdta[i]<<endl;
	}
	if(!author.empty())
	{
		for(i=0;i<(int)author.size();i++)
			file<<author[i]<<endl;
	}
	if(!revdat.empty())
	{
		for(i=0;i<(int)revdat.size();i++)
			file<<revdat[i]<<endl;
	}
	if(!sprsde.empty())
	{
		for(i=0;i<(int)sprsde.size();i++)
			file<<sprsde[i]<<endl;
	}
	if(!jnrl.empty())
	{
		for(i=0;i<(int)jnrl.size();i++)
			file<<jnrl[i]<<endl;
	}
	if(!jrnl.empty())
	{
		for(i=0;i<(int)jrnl.size();i++)
			file<<jrnl[i]<<endl;
	}
	if(!remark.empty())
	{
		for(i=0;i<(int)remark.size();i++)
			file<<remark[i]<<endl;
	}
	if(!dbref.empty())
	{
		for(i=0;i<(int)dbref.size();i++)
			file<<dbref[i]<<endl;
	}
	if(!seqadv.empty())
	{
		for(i=0;i<(int)seqadv.size();i++)
			file<<seqadv[i]<<endl;
	}
	if(!seqres.empty())
	{
		for(i=0;i<(int)seqres.size();i++)
			file<<seqres[i]<<endl;
	}
	if(!modres.empty())
	{
		for(i=0;i<(int)modres.size();i++)
			file<<modres[i]<<endl;
	}
	if(!ftnote.empty())
	{
		for(i=0;i<(int)ftnote.size();i++)
			file<<ftnote[i]<<endl;
	}
	if(!het.empty())
	{
		for(i=0;i<(int)het.size();i++)
			file<<het[i]<<endl;
	}
	if(!hetnam.empty())
	{
		for(i=0;i<(int)hetnam.size();i++)
			file<<hetnam[i]<<endl;
	}
	if(!hetsyn.empty())
	{
		for(i=0;i<(int)hetsyn.size();i++)
			file<<hetsyn[i]<<endl;
	}
	if(!formul.empty())
	{
		for(i=0;i<(int)formul.size();i++)
			file<<formul[i]<<endl;
	}
	if(!helix.empty())
	{
		for(i=0;i<(int)helix.size();i++)
			file<<helix[i]<<endl;
	}
	if(!sheet.empty())
	{
		for(i=0;i<(int)sheet.size();i++)
			file<<sheet[i]<<endl;
	}
	if(!turn.empty())
	{
		for(i=0;i<(int)turn.size();i++)
			file<<turn[i]<<endl;
	}
	if(!ssbond.empty())
	{
		for(i=0;i<(int)ssbond.size();i++)
			file<<ssbond[i]<<endl;
	}
	if(!hydbnd.empty())
	{
		for(i=0;i<(int)hydbnd.size();i++)
			file<<hydbnd[i]<<endl;
	}
	if(!sltbrg.empty())
		{
		for(i=0;i<(int)sltbrg.size();i++)
			file<<sltbrg[i]<<endl;
	}	
	if(!link.empty())
	{
		for(i=0;i<(int)link.size();i++)
			file<<link[i]<<endl;
	}
	if(!cispep.empty())
	{
		for(i=0;i<(int)cispep.size();i++)
			file<<cispep[i]<<endl;
	}
	if(!site.empty())
	{
		for(i=0;i<(int)site.size();i++)
			file<<site[i]<<endl;
	}
	if(!cryst1.empty())
	{
		for(i=0;i<(int)cryst1.size();i++)
			file<<cryst1[i]<<endl;
	}
	if(!origx1.empty())
	{
		for(i=0;i<(int)origx1.size();i++)
			file<<origx1[i]<<endl;
	}
	if(!origx2.empty())
	{
		for(i=0;i<(int)origx2.size();i++)
			file<<origx2[i]<<endl;
	}
	if(!origx3.empty())
	{
		for(i=0;i<(int)origx3.size();i++)
			file<<origx3[i]<<endl;
	}
	if(!scale1.empty())
	{
		for(i=0;i<(int)scale1.size();i++)
			file<<scale1[i]<<endl;
	}
	if(!scale2.empty())
	{
		for(i=0;i<(int)scale2.size();i++)
			file<<scale2[i]<<endl;
	}
	if(!scale3.empty())
	{
		for(i=0;i<(int)scale3.size();i++)
			file<<scale3[i]<<endl;
	}
	if(!mtrix1.empty())
	{
		for(i=0;i<(int)mtrix1.size();i++)
			file<<mtrix1[i]<<endl;
	}
	if(!mtrix2.empty())
	{
		for(i=0;i<(int)mtrix2.size();i++)
			file<<mtrix2[i]<<endl;
	}
	if(!mtrix3.empty())
	{
		for(i=0;i<(int)mtrix3.size();i++)
			file<<mtrix3[i]<<endl;
	}
	if(!tvect.empty())
	{
		for(i=0;i<(int)tvect.size();i++)
			file<<tvect[i]<<endl;
	}
}
//anisou��ʱδ����
void Protein::OutputSucField(ofstream& file)
{
	int i;
	if(!hetatm.empty())
	{
			for(i=0;i<(int)hetatm.size();i++)
				file<<hetatm[i]<<endl;
	}
	if(!conect.empty())
	{
			for(i=0;i<(int)conect.size();i++)
				file<<conect[i]<<endl;
	}
	if(!master.empty())
	{
			for(i=0;i<(int)master.size();i++)
				file<<master[i]<<endl;
	}
	if(!end.empty())
	{
		for(i=0;i<(int)end.size();i++)
				file<<end[i]<<endl;
	}
	else
		file<<"END                                                                             ";

			
}
/**************************************************
 * ��ÿ��������Ϣ�ֱ�����ļ���
 * ���filenameΪĬ��ֵ�����ļ���Ϊpdb��ʶ+����ʶ.ent
 * ����ļ�����ΪĬ��ֵ��������ʶ���뵽���һ��.֮ǰ
 * ����ǰ�ĵ������н���һ���������ļ����в�������ʶ
 * �������ʶΪ�ո�������Ҫ��ʾʱ����A��ʾ
 ***************************************************/
bool Protein::SaveChains(string dirname,string filename)
{
	
	string line,fname,tempstr;
	Atom atom;
	int j,po;
	int num_atom;
	int chain_num;
	if(num_chain<1)
		return false;

	for(chain_num=0;chain_num<num_chain;chain_num++)
	{
		//�����ļ���
		if(filename!="AUTO")
		{
			fname=filename;
			if(chain_num>1)
			{
				if(m_chains[chain_num].chain_identifier!=' ')
				{
					po=(int)fname.find('.');
					tempstr=m_chains[chain_num].chain_identifier;
					fname.insert(po,tempstr);
				}
				else
				{
					po=(int)fname.find('.');
					fname.insert(po,"A");
				}
			}
		}
		else
		{
			fname=PDB_id;
			if(num_chain>1)
			{
				if(m_chains[chain_num].chain_identifier!=' ')
				{					
					fname+=m_chains[chain_num].chain_identifier;
				}
				else
					fname+='A';
			}
			fname+=".ent";
		}
		//���Ŀ¼��ʶ
		if(dirname!="AUTO")
		{
			_mkdir(dirname.c_str());
			fname=dirname+"\\"+fname;
		}
		ofstream file;		
		file.open(fname.c_str(),ios::out);
		if(!file.is_open())
			return false;
		/*line="HEADER";
		line.resize(80);
		for(i=0;i<(int)PDB_id.length();i++)
			line[62+i]=PDB_id[i];
		file<<line<<endl;*/
		
		OutputPreField(file);		
		num_atom=m_chains[chain_num].num_atom;
		for(j=0;j<num_atom;j++)
		{
			atom=m_chains[chain_num].m_atoms[j];
			Atom2Line(atom,line);
			file<<line<<endl;
		}
		Atom2Ter(atom,line);
		file<<line<<endl;
		OutputSucField(file);
	}
	return true;
}
/**********************************************
 * ͳ�ư�����ԹǼ�Ťת�ǵ�������
 * �Ǽ�Ťת�Ǳ����Լ���
 * ����Ǽ�Ťת�Ƿֳ�36*36��
 * �ɹ�ʱ����true�����򷵻�false
 *********************************************/
bool Chain::Pre_Dihedral(vector<vector<double> >& preference)
{
	int len,i,bin,bpsi,bphi,aa;
	double phi,psi;
	len=(int)m_residues.size();
	for(i=0;i<len;i++)
	{
		//����Ǽ�Ťת������
		phi=m_residues[i].phi;
		psi=m_residues[i].psi;
		if(fabs(phi+999)<1e-5)//����ǲ�������ͳ��
			continue;
		if(fabs(psi+999)<1e-5)
			continue;
		bphi=(int)((phi+180)/10);
		if(bphi>=36)
			bphi=36;
		bpsi=(int)((psi+180)/10);
		if(bpsi>=36)
			bpsi=36;
		bin=bphi*36+bpsi;
		if(bin>1295)
			return false;
		//���㰱�������
		aa=Rname2ID(m_residues[i].c_name);		
		preference[bin][aa]+=1.0;
	}
	return true;
}
//����ÿ��������ĽӴ�����
bool Chain::Cal_Num_Contact()
{
	int i,j,len;
	int index,atom_index;
	Point pi,pj;
	double distance;
	len=num_residue;

	for(i=0;i<len;i++)
	{
		index=Atomname2ID("CA");
		atom_index=m_residues[i].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"����Ӵ�����ʱ�����ڰ�����"<<i+1<<"���Ҳ���CAԭ��"<<endl;
			return false;
		}
		pi=m_atoms[atom_index].pt;
		for(j=i+1;j<len;j++)
		{
			index=Atomname2ID("CA");
			atom_index=m_residues[j].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"����Ӵ�����ʱ�����ڰ�����"<<j+1<<"���Ҳ���CAԭ��"<<endl;
				return false;
			}
			pj=m_atoms[atom_index].pt;
			distance=Distance(pi,pj);
			if(distance<=8)//����С��8��������Ϊ����contact
			{
				m_residues[i].num_contact++;
				m_residues[j].num_contact++;
			}
		}
	}
	return true;
}
//ͳ�ư�����ԽӴ��������͵ĳ��ִ������Ӵ����������Ѿ�����
bool Chain::Pre_Num_Contact(vector<vector<double> >& preference)
{
	int i,len;
	int bin,aa;
	len=num_residue;

	for(i=0;i<len;i++)
	{
		bin=m_residues[i].num_contact;
		if(bin<1)
			bin=1;
		else if(bin>25)
			bin=25;
		bin-=1;//ת������
		aa=Rname2ID(m_residues[i].c_name);	
		preference[bin][aa]+=1.0;
	}
	return true;
}
bool IsDigit(string str)
{
	int i,len;
	len=(int)str.size();
	for(i=0;i<len;i++)
	{
		if((!isdigit(str[i]))&&(str[i]!=' ')&&(str[i]!='-'))
			return false;
	}
	return true;
}
/*********************************************************************
 * �ۼ�ԭ������1������
 * ԭ������1������ָ��ԭ���������ϵļ������k�����У���11��22Ϊ�ֽ��
 * ԭ���ڿռ�ľ�������0�仯��10����0.5Ϊ�����20��
 * �����������ϵ���ԭ�Ӷ�25��
 * ��Ӧ�����������ĳ���Ϊ1500
 *********************************************************************/
bool Chain::Pre_Num_Atomtype1(vector<double>& preference)
{
	int i1,i2,len;
	int k,s,daa,po,kk,type1,type2;
	Point p1,p2;
	double distance;
	string atomname1,atomname2;

	len=this->num_atom;
	for(i1=0;i1<len;i1++)
	{
		atomname1=m_atoms[i1].name;
		if((atomname1=="CB")||(atomname1=="N")||(atomname1=="O")||(atomname1=="CA")||(atomname1=="C"))//ԭ�������ڶ���ʱ�Ѿ�ת��Ϊ��д
		{
			p1=m_atoms[i1].pt;
			for(i2=i1+1;i2<len;i2++)
			{
				atomname2=m_atoms[i2].name;
				if((atomname2=="CB")||(atomname2=="N")||(atomname2=="O")||(atomname2=="CA")||(atomname2=="C"))
				{
					p2=m_atoms[i2].pt;
					distance=Distance(p1,p2);
					if(distance>=10.0)//�����þ����ԭ�ӶԺ��Ե�
						continue;
					//ȷ���ռ������s
					s=int(distance*2);
					if(s>=20)
						s=19;
					//ȷ�������������k
					kk=m_atoms[i1].residule_chain_seqm-m_atoms[i2].residule_chain_seqm;
					if(kk<0)
						kk=-kk;
					if(kk<11)
						k=0;
					else if(kk>=11||kk<=22)
						k=1;
					else
						k=2;
					//ȷ��ԭ�ӶԱ��
					if(atomname1=="CA")
						type1=0;
					else if(atomname1=="CB")
						type1=1;
					else if(atomname1=="N")
						type1=2;
					else if(atomname1=="C")
						type1=3;
					else //"O"
						type1=4;
					if(atomname2=="CA")
						type2=0;
					else if(atomname2=="CB")
						type2=1;
					else if(atomname2=="N")
						type2=2;
					else if(atomname2=="C")
						type2=3;
					else //"O"
						type2=4;
					daa=type1*5+type2;
					po=k*500+s*25+daa;
					preference[po]+=1.0;
				}
			}
		}
	}
	return true;
}
/*********************************************************************
 * �ۼ�ԭ������2������
 * ԭ������1������ָ��ԭ���������ϵļ������k�����У���11��22Ϊ�ֽ��
 * ԭ���ڿռ�ľ�������1.5�仯��11.5����1Ϊ�����10��
 * ԭ�����͹�40��,ԭ�Ӷ�1600��
 * ��Ӧ�����������ĳ���Ϊ48000
 *********************************************************************/
bool Chain::Pre_Num_Atomtype2(vector<double>& preference)
{
	int i1,i2,len;
	int k,s,daa,po,kk,type1,type2;
	Point p1,p2;
	double distance;
	string atomname1,atomname2;	
	char c1,c2;
	int aaindex,atomindex;

	len=this->num_atom;
	for(i1=0;i1<len;i1++)
	{
		atomname1=m_atoms[i1].name;
		c1=m_atoms[i1].c_resname;
		aaindex=Rname2ID20(c1);//��ԭ�����ڵİ���������ת����20����׼���������
		if(aaindex==-1)
		{
			cout<<"����:��������Ч�İ�������ŵ�"<<i1<<"��ԭ�����ڵİ���������Ϊ"<<c1<<endl;
			cout<<"��ԭ�ӽ�������"<<endl;
			continue;
		}
		atomindex=Atomname2ID(atomname1.c_str());//ԭ������ת�������
		if(atomindex==36)//����OXTԭ��
			continue;
		if(atomindex==-1)
		{
			cout<<"����:��������Ч��ԭ����ŵ�"<<i1<<"��ԭ��������Ϊ"<<atomname1.c_str()<<endl;
			cout<<"��ԭ�ӽ�������"<<endl;
			continue;
		}
		type1=Atomtype40[aaindex][atomindex];
		if(type1==-1)
		{
			cout<<"����:δ�ܻ��ԭ�ӵ���ȷ��Ű���������"<<c1<<" ���"<<aaindex<<" ԭ������"<<atomname1.c_str()<<" ���"<<atomindex<<endl;
			cout<<"��ԭ�ӽ�������"<<endl;
			continue;
		}
		//if((atomname1=="CB")||(atomname1=="N")||(atomname1=="O")||(atomname1=="CA")||(atomname1=="C"))//ԭ�������ڶ���ʱ�Ѿ�ת��Ϊ��д
		{
			p1=m_atoms[i1].pt;
			for(i2=i1+1;i2<len;i2++)
			{
				atomname2=m_atoms[i2].name;
				c2=m_atoms[i2].c_resname;
				aaindex=Rname2ID20(c2);//��ԭ�����ڵİ���������ת����20����׼���������
				if(aaindex==-1)
				{
					cout<<"����:��������Ч�İ�������ŵ�"<<i2<<"��ԭ�����ڵİ���������Ϊ"<<c2<<endl;
					cout<<"��ԭ�ӽ�������"<<endl;
					continue;
				}
				atomindex=Atomname2ID(atomname2.c_str());//ԭ������ת�������
				if(atomindex==36)//����OXTԭ��
					continue;
				if(atomindex==-1)
				{
					cout<<"����:��������Ч��ԭ����ŵ�"<<i2<<"��ԭ��������Ϊ"<<atomname2.c_str()<<endl;
					cout<<"��ԭ�ӽ�������"<<endl;
					continue;
				}
				type2=Atomtype40[aaindex][atomindex];
				if(type2==-1)
				{
					cout<<"����:δ�ܻ��ԭ�ӵ���ȷ��Ű���������"<<c2<<" ���"<<aaindex<<" ԭ������"<<atomname2.c_str()<<" ���"<<atomindex<<endl;
					cout<<"��ԭ�ӽ�������"<<endl;
					continue;
				}
				//if((atomname2=="CB")||(atomname2=="N")||(atomname2=="O")||(atomname2=="CA")||(atomname2=="C"))
				{
					p2=m_atoms[i2].pt;
					distance=Distance(p1,p2);
					if((distance<1.5)||(distance>=11.5))//�����þ����ԭ�ӶԺ��Ե�
						continue;
					//ȷ���ռ������s
					s=(int)(distance-1.5);
					if(s>=10)
						s=9;
					//ȷ�������������k
					kk=m_atoms[i1].residule_chain_seqm-m_atoms[i2].residule_chain_seqm;					
					if(kk<0)
						kk=-kk;
					if(kk<=1) //���ڻ���ͬһ��������֮�ڵĲ�����
						continue;
					if(kk<11)
						k=0;
					else if(kk>=11||kk<=22)
						k=1;
					else
						k=2;					
					daa=(type1-1)*40+(type2-1);
					po=k*16000+s*1600+daa;
					preference[po]+=1.0;
				}
			}
		}
	}
	return true;
}






