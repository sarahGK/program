#include "stdafx.h"
#include "hr_Info.h"
#include "rd_pdb.h"
#include "rmsd.h"


HR_genel::HR_genel(){
	num_decoys = 500;
	
}

HR_genel::~HR_genel(){}

void HR_genel::set_prbrank(){
	vector<int> decoystemp;
	for(int i=0;i<num_decoys;i++)
		decoystemp.push_back(decoys[i]);
	for(int i=0;i<num_decoys-1;i++){
		int p = i;
		int t;
		double d;
		for(int j=i+1;j<num_decoys;j++)
			if(probal2_decoys[j]>probal2_decoys[p]) p=j;
		prbrank2.push_back(decoystemp[p]);
		if(p!=i){
			d = probal2_decoys[p];
			probal2_decoys[p] = probal2_decoys[i];
			probal2_decoys[i] = d;
			t = decoystemp[p];
			decoystemp[p] = decoystemp[i];
			decoystemp[i] = t;
		}
	}
		prbrank2.push_back(decoystemp[num_decoys-1]);
}

double HR_genel::cpt_rcgrate(int no){
	if(no==1){
		if(prbvalue1_native>probal1_decoys[0]) recognition_rate= 1;
		else recognition_rate= 0;
	}
	else{
		if(prbvalue2_native>probal2_decoys[0]) recognition_rate= 1;
		else recognition_rate= 0;
	}
	return recognition_rate;
}

int HR_genel::get_rank(int no){
	int i=0;
	if(no==1)	while(i<num_decoys&&prbvalue1_native<probal1_decoys[i])		i++;
	else  while(i<num_decoys&&prbvalue2_native<probal2_decoys[i])		i++;
	rank = i+1;
	return rank;
}

double HR_genel::get_z(int no){
	double sum1=0.0;
	double sum2=0.0;
	double denominator ;
	
	for(int i=0;i<num_decoys;i++){
		if(no==1){
		sum1 = sum1 + probal1_decoys[i];
		sum2 = sum2 + probal1_decoys[i]*probal1_decoys[i];
		}
		else{
		sum1 = sum1 + probal2_decoys[i];
		sum2 = sum2 + probal2_decoys[i]*probal2_decoys[i];
		}
	}
	sum1 = sum1/num_decoys;
	sum2 = sum2/num_decoys;
	denominator = sqrt(sum2-sum1*sum1);
	if(no==1) z_score = (sum1- prbvalue1_native)/denominator;
	else z_score = (sum1-prbvalue2_native)/denominator;
	return z_score;
}

double HR_genel::get_z(){
	double sum1=0.0;
	double sum2=0.0;
	double denominator ;
	for(int i=0;i<num_decoys;i++){	
			sum1 = sum1 + probal2_decoys[i];
			sum2 = sum2 + probal2_decoys[i]*probal2_decoys[i];
	}
	sum1 = sum1/num_decoys;
	sum2 = sum2/num_decoys;
	denominator = sqrt(sum2-sum1*sum1);
     z_score = (sum1-prbvalue2_native)/denominator;
	return z_score;
}

int HR_genel::find_mark(int x){
	for (int i=0;i<num_decoys;i++)
		if (x==decoys[i]) return i;
}

double HR_genel::get_rel(int no){
	double nominator = 0.0;
	double denominator = 0.0;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double ave1, ave2;
	for(int i=0;i<num_decoys;i++){
		if(no==1)		sum1 = sum1 + probal1_decoys[i];
		else 	sum1 = sum1 + probal2_decoys[i];
		sum2 = sum2 + rmsd[i];
	}
	ave1 = sum1/num_decoys;
	ave2 = sum2/num_decoys;
	sum1 = sum2 = 0.0;
	int index;
	for(int i=0;i<num_decoys;i++){
		if(no==1) {
			index = find_mark(prbrank1[i]);
			nominator = nominator + (probal1_decoys[i]-ave1)*(rmsd[index]-ave2);
			sum1 = sum1+ (probal1_decoys[i]-ave1) * (probal1_decoys[i]-ave1);
		}
		else{
			index = find_mark(prbrank2[i]);
			nominator = nominator + (probal2_decoys[i]-ave1)*(rmsd[index]-ave2);
			sum1 = sum1+ (probal2_decoys[i]-ave1) * (probal2_decoys[i]-ave1);
		}
		sum2 = sum2 + (rmsd[i]-ave2)* (rmsd[i]-ave2);
		
	}
	denominator = sqrt(sum1 * sum2);
	revalance= nominator/denominator;
	return revalance;
}

double HR_genel::get_rel(){
	double nominator = 0.0;
	double denominator = 0.0;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double ave1, ave2;
	for(int i=0;i<num_decoys;i++){
		sum1 = sum1 + probal2_decoys[i];
		sum2 = sum2 + rmsd[i];
	}
	ave1 = sum1/num_decoys;
	ave2 = sum2/num_decoys;
	sum1 = sum2 = 0.0;
	int index;
	for(int i=0;i<num_decoys;i++){	
			index = find_mark(prbrank2[i]);
			nominator = nominator + (probal2_decoys[i]-ave1)*(rmsd[index]-ave2);
			sum1 = sum1+ (probal2_decoys[i]-ave1) * (probal2_decoys[i]-ave1);
		    sum2 = sum2 + (rmsd[i]-ave2)* (rmsd[i]-ave2);
	}
	denominator = sqrt(sum1 * sum2);
	revalance= nominator/denominator;
	return revalance;
}

double HR_genel::count_rmsd(int no){
	double result = 0.0;
	string filename1,filename2;
	filename1 = "G:\\hr\\hrdecoys_selected\\";
	filename1.append(pdbid.c_str(),5);
	filename1.append("\\",1);
	filename1.append(pdbid.c_str(),5);
	filename1.append(".",1);
	filename2 = filename1;
	char ch[5]="" ;
	if (no==1)  filename1.append(itoa(prbrank1[0],ch,10),4) ;
	else filename1.append(itoa(prbrank2[0],ch,10),4) ;
	filename1=filename1.c_str();
	filename1.append(".coord.pdb",10);
	double min  = rmsd[0];
	int index = 0;
	for (int i=1;i<num_decoys;i++)
		if (rmsd[i]<min){
			min = rmsd[i];
			index = i;
		}
	filename2.append(itoa(decoys[index],ch,10),4) ;
	filename2 = filename2.c_str();
	filename2.append(".coord.pdb",10);
	if(Ev_Rmsd(filename1,filename2,1,result))
		return result;
	return 1234567;
}

int HR_genel::find(int x){
	for (int i=0;i<N*num_decoys;i++)
	{
		if (x==decoys[i]) return 1;
	}
	return 0;
}

double HR_genel::get_nenrich(int no){
	double d;
	int temp;
	double sum = 0.0;
	for (int i =0;i<N*num_decoys;i++)
	{
		int p = i;
		for(int j=i+1;j<num_decoys;j++)
			if (rmsd[j]<rmsd[p]) p = j;
		if (p!=i)
		{
			d = rmsd[p];
			rmsd[p] = rmsd[i];
			rmsd[i] = d;
			temp = decoys[p];
			decoys[p] = decoys[i];
			decoys[i] = temp;
		}
	}

	for (int i=0;i<N*num_decoys;i++){
			if (no==1)  sum = sum + find(prbrank1[i]);
			else sum = sum + find(prbrank2[i]);
		}
	sum = sum/(N*num_decoys)/N;
	return sum;
}

double HR_genel::get_re(int no){
	double nominator = 0.0;
	double denominator = 0.0;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double ave1, ave2;
	for(int i=0;i<num_decoys;i++){
		if(no==1)		sum1 = sum1 + probal1_decoys[i];
		else 	sum1 = sum1 + probal2_decoys[i];
		sum2 = sum2 + rmsd[i];
	}
	ave1 = sum1/num_decoys;
	ave2 = sum2/num_decoys;
	sum1 = sum2 = 0.0;

	for(int i=0;i<num_decoys;i++){
		if(no==1) {
			nominator = nominator + (probal1_decoys[i]-ave1)*(rmsd[i]-ave2);
			sum1 = sum1+ (probal1_decoys[i]-ave1) * (probal1_decoys[i]-ave1);
		}
		else{
			
			nominator = nominator + (probal2_decoys[i]-ave1)*(rmsd[i]-ave2);
			sum1 = sum1+ (probal2_decoys[i]-ave1) * (probal2_decoys[i]-ave1);
		}
		sum2 = sum2 + (rmsd[i]-ave2)* (rmsd[i]-ave2);

	}
	denominator = sqrt(sum1 * sum2);
	revalance= nominator/denominator;
	return revalance;
}
