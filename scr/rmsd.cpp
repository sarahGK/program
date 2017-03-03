#include "stdafx.h"
#include "rmsd.h"

/***************************************************************
 * 评估两个蛋白质的RMSD,仅对第一个链比较
 * 成功时返回true
 * 否则返回false,此时,可能是文件不存在,或者序列不一致
 * leibie指定评估的类型:1:仅Ca原子2:所有原子其他返回false
 * rmsd接受结果
 **************************************************************/
bool Ev_Rmsd(const string& filename1,const string& filename2, int leibie,double& rmsd)
{
	int i;
	Protein p1,p2;
	vector<Point> vp1,vp2;
	if((leibie!=1)&&(leibie!=2))
		return false;
	if(!p1.LoadFromPDB(filename1))
		return false;
	if(!p2.LoadFromPDB(filename2))
		return false;
	if(p1.num_chain<=0)
		return false;
	if(p2.num_chain<=0)
		return false;
	if(p1.m_chains[0].seq!=p2.m_chains[0].seq)
		return false;
	int len;
	if(leibie==1)
		len=(int)p1.m_chains[0].m_residues.size();
	else
		len=(int)p1.m_chains[0].m_atoms.size();
	vp1.resize(len);
	vp2.resize(len);
	int index,atom_index;
	for(i=0;i<len;i++)
	{
		if(leibie==1)
		{
			index=Atomname2ID("CA");
			atom_index=p1.m_chains[0].m_residues[i].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"计算RMSD时出错：在蛋白质的p1氨基酸"<<i<<"中找不到CA原子"<<endl;
				return false;
			}
			vp1[i]=p1.m_chains[0].m_atoms[atom_index].pt;

			index=Atomname2ID("CA");
			atom_index=p2.m_chains[0].m_residues[i].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"计算RMSD时出错：在蛋白质p2的氨基酸"<<i<<"中找不到CA原子"<<endl;
				return false;
			}
			vp2[i]=p2.m_chains[0].m_atoms[atom_index].pt;
		}
		else
		{
			vp1[i]=p1.m_chains[0].m_atoms[i].pt;
			vp2[i]=p2.m_chains[0].m_atoms[i].pt;
		}
	}
	double u[3][3],t[3];
	double result;
	result=kabsch(vp1,vp2,u,t);
	if(result<0)
		return false;
	rmsd=result;
	return true;
}

	

