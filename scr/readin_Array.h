//const char rname_one_rank20[21] = "ACDEFGHIKLMNPQRSTVWY";
double proposal1[20][1296];//task2�е�һ�ַ���������ĸ��ʱ�ĸ���ֵ
double proposal2[20][15];//task2�е�һ�ַ���������ĸ��ʱ�ĸ���ֵ

/*int BinarySearchChar(char c)
{
	int left,right,middle;
	left=0;
	right=20;
	if(c<rname_one_rank20[0])
		return -1;
	if(c>rname_one_rank20[19])
		return -1;
	while(right-left>1)
	{
		middle=(left+right)/2;
		if(c==rname_one_rank20[middle])
			return middle;
		if(c<rname_one_rank20[middle])
			right=middle;
		else
			left=middle;
	}
	if(c==rname_one_rank20[left])
		return left;
	else if(c==rname_one_rank20[right])
		return right;
	else 
		return -1;
}*/



void initial(){//����task2�ĸ�������
	ifstream datain;
	datain.open("E:\\talk\\biology\\task2\\probalistics.txt",ios::in);
	//datain.open("G:\\task3\\result\\probability2.txt",ios::in);
	if(!datain.is_open()){
		cout<<"can't open the file!"<<endl;
		return;
	}
	cout<<"���ڶ�ȡ�ļ���"<<"E:\\talk\\biology\\task2\\probalistics.txt"<<endl;
	//cout<<"���ڶ�ȡ�ļ���"<<"G:\\task3\\result\\probability2.txt"<<endl;
	for(int i=0;i<20;i++)
		for(int j=0;j<1296;j++)
			datain>>proposal1[i][j];
	for(int i=0;i<20;i++)
		for(int j =0;j<15;j++)
			datain>>proposal2[i][j];
}
/*
*�Զ���ķ�ʽ��ָ�����ļ�
*/
bool openfile(string filename,ifstream & datain){
	
	datain.open (filename.c_str(),ios::in);
	if(!datain.is_open()){
		printf("�޷����ļ���%s\n",filename.c_str());
		return false;
	}
	else{
		cout<<"���ڶ�ȡ�ļ���"<<filename.c_str()<<endl;
		return true;
	}
}
/*
*task2�жԷ���һÿ����ĸ���ֵĸ���ֵȡ����
*/
double log_probal1(char ch,int x){
	
	int i = BinarySearchChar(ch);
	if(i!=-1&&x>=0&&x<1296)	return log(proposal1[i][x]);
	else return 0;
}
/*
*task2�жԷ�����ÿ����ĸ���ֵĸ���ֵȡ����
*/
double log_probal2(char ch,int x){
	int i = BinarySearchChar(ch);
	if(i!=-1&&x>=1&&x<=15)	{
		//if(proposal2[i][x-1]<=0)
		//	double p = proposal2[i][x-1];
		return log(proposal2[i][x-1]);
	}
	else return 0;
}
/*
*task3�жԷ�����ÿ����ĸ���ֵĸ���ֵȡ����
*/
double log_probal2(char c1,int x1,char c2,int x2,counter & c){
	int i = BinarySearchChar(c1);
	int j = BinarySearchChar(c2);
	if(i!=-1&&j!=-1&&x1>=1&&x1<=15&&x2>=1&&x2<=15) {
		//if(c.prob[i][x1-1][j][x2-1]<=0)
		//	double p = c.aprob[i][x1-1][j][x2-1];
		return log(c.aprob[i][x1-1][j][x2-1]);
	}
	else return 0;
}

/*void print_Arr(int length){
	for(int i=0;i<20;i++){
		for(int j=0;j<length;j++){
			if(length>15)	cout<<proposal1[i][j]<<"  ";
			else	cout<<proposal2[i][j]<<"  ";
		}
		cout<<endl;		
	}
}*/
/*
*��hr���ļ�����ĸ����Ϣ����allhr��
*/
void prod_hrall(vector<HR_genel> & allhr,counter & c){
    
	ifstream name;
	string nameline;
	int i = 1;
	
	if(openfile("E:\\talk\\biology\\task1\\alldecoys.txt",name)){//���ռ�¼�Ŷ��ṹ�ļ�����˳��
			HR_genel hr;
			bool flag1=false;//�Ƿ������Ȼ�ṹ
			bool flag2=false;//�Ƿ����500���Ŷ��ṹ
			bool flag3=false;//�Ƿ�����Ŷ��ṹ�ı�ʾ
		while(!name.eof()){
		
			if(i==1){//�����ļ���һ��
				i++;
				getline(name,nameline);			
				continue;
			}
			if(i%2==0){//ȷ����ȡ����ǵ��ļ���
				string dihedralfile = "G:\\hr\\result\\";
				getline(name,nameline);				
				dihedralfile.append(nameline.substr(0,5).c_str(),5);
				dihedralfile.append(".txt",4);
				ifstream dihedraldata;
				if(openfile(dihedralfile,dihedraldata)){
					hr.set_ID(nameline.substr (0,5));
					string dihedralline;
					getline(dihedraldata,dihedralline);
				
					
					while(!dihedraldata.eof()){
				
					if(dihedralline.length()==17){//��ǰ��ȡ������Ȼ�ṹ����Ϣ
						flag1 = true;
						cout<<"���ڼ�����Ȼ�ṹ�ĸ���ֵ������"<<endl;
						getline(dihedraldata,dihedralline);
						double sum1,sum2;
						sum1= sum2 = 0.0;
						bool first = true;//��ʾ�Ƿ��ǵ�һ�����������
						int x;//ǰ��������������ĸ
						char c1;//��¼ǰ���������ĸ��ʾ
						while(dihedralline.length()>20){
						char ch = dihedralline.substr(8,1).c_str()[0];
						//int x = atoi(dihedralline.substr(30,10).c_str());
						int y = atoi(dihedralline.substr(40,12).c_str());
						//sum1 = sum1 + log_probal1(ch,x);
						if(first){
							first = false;
							sum2 = sum2 + log_probal2(ch,y);//��һ����������ӣ�����task2�������ֵ
						}else sum2 = sum2 + log_probal2(c1,x,ch,y,c);//���ఱ���ᰴ��task3�������ֵ
						c1 = ch;//��ǰֵΪ��һ�������ǰ���������ĸ��ʾ
						x = y;//��ǰֵΪ��һ�������ǰ����������ǵ���ĸ��ʾ
						getline(dihedraldata,dihedralline);
						}
						//hr.set_prbvalue(sum1,sum2);
						hr.set_prbvalue(sum2);
			
					}
					
					if(flag1){//��ȡ����Ȼ�ṹ����Ϣ�󣬶�ȡ�Ŷ��ṹ����Ϣ
					
						cout<<"���ڼ����Ŷ��ṹ�ĸ���ֵ������"<<endl;
						getline(dihedraldata,dihedralline);
						double sum1 ,sum2;
						sum1= sum2 = 0.0;
						bool first = true;
						int x;
						char c1;
						while(!dihedraldata.eof()&&dihedralline.length()>21){
						
						char ch = dihedralline.substr(8,1).c_str()[0];
						//int x = atoi(dihedralline.substr(30,10).c_str());
						int y = atoi(dihedralline.substr(40,12).c_str());
			
						//sum1 = sum1 + log_probal1(ch,x);
						if (first)
						{
							first = false;
							sum2 = sum2 + log_probal2(ch,y);
						}else sum2 = sum2 + log_probal2(c1,x,ch,y,c);
						x = y;
						c1 = ch;
						getline(dihedraldata,dihedralline);
						}
						//hr.set_pdecoys(sum1,sum2);		
						hr.set_pdecoys(sum2);
					}
					flag2 = true;
					}				
				}
					i++;
					continue;
			}
		
			else{//alldecoys.txt�������м�¼��ȡ�Ŷ��ṹ�ı�ʾ
				cout<<"���ڶ�ȡ�Ŷ��ṹ�ı�ʾ������"<<endl;
				int m[500] ;
				for(int k=0;k<500;k++){	
				name>>m[k];
				hr.set_mark(m[k]);
				}
				getline(name,nameline);
				i++;
				flag3 = true;
			}
			
			if(flag1&&flag2&&flag3)	{//��ÿһ�������ʵ���Ȼ�ṹ���Ŷ��ṹ�Լ���ʾ����ȡ��ϣ�����ǰ��hr��ӵ�������
				allhr.push_back(hr);	
				hr.~HR_genel();
			}
			cout<<allhr.size()<<endl;
			if(!flag3||!flag1||!flag2) continue;
			//HR_genel hr;
			flag1=flag2=flag3=false;//������һ�������ʷ���
		}
	}
}

/*
*�����Ŷ��ṹ��rmsdֵ
*/
void rmsd_allhr(vector<HR_genel> & allhr){
	for(int i=0;i<allhr.size();i++){
		allhr[i].set_prbrank();//��ÿһ��hr�ĸ���ֵ����
		 ifstream rmsdfile;
		 string filename="G:\\hr\\hrdecoys_selected\\";
		 filename.append(allhr[i].pdbid.c_str(),5);
		 filename.append("\\",1);
		 filename.append(allhr[i].pdbid.c_str(),5);
		 filename.append(".rmsd_recal",11);
		 if(openfile(filename,rmsdfile)){
				cout<<"���ڶ����Ŷ��ṹ��rmsdֵ������"<<endl;
				string line;
				int m;
				double d;
				for(int j=0;j<allhr[i].decoys.size();j++){
					rmsdfile>>m;					
					while(!rmsdfile.eof()&&allhr[i].decoys[j]!=m){
						getline(rmsdfile,line);
						rmsdfile>>m;
					}
					rmsdfile>>d;
					allhr[i].set_rmsd(d);
					getline(rmsdfile,line);
				}
		 }
	}
}
/*
*����6������ָ��ֵ�����е����ʵ�ƽ��ֵ,�������ָ���ļ�
*/
void count6(vector<HR_genel> & allhr,ofstream & of1,ofstream & of2){
	//double sumofrate1=0.0;
	double sumofrate2=0.0;
	//double sumofrank1=0.0;
	double sumofrank2=0.0;
	//double sumofrlv1 = 0.0;
	double sumofrlv2 = 0.0;
	//double sumofzs1 = 0.0;
	double sumofzs2 = 0.0;
	//double sumofdrmsd1 = 0.0;
	double sumofdrmsd2 = 0.0;
	//double sumofenrich1 = 0.0;
	double sumofenrich2 = 0.0;
	of1<<"pdbid"<<"   "<<"ʶ����"<<"   "<<"��"<<"   "<<"�����"<<"   "<<"z����"<<"   "<<"drmsd"<<"   "<<"Nenrichment"<<endl;
	of2<<"pdbid"<<"   "<<"ʶ����"<<"   "<<"��"<<"   "<<"�����"<<"   "<<"z����"<<"   "<<"drmsd"<<"   "<<"Nenrichment"<<endl;
	of1.width(9);
	of2.width(9);
	for(int i=0;i<allhr.size();i++){
		of1<<allhr[i].pdbid<<"   ";
		of2<<allhr[i].pdbid<<"   ";
		double recgnition = 0.0;
		/*recgnition = allhr[i].cpt_rcgrate(1);
		sumofrate1 = sumofrate1 + recgnition;
		cout<<"��"<<i<<"�������ʽṹ��ʶ����Ϊ��"<<recgnition<<"  ";
		of1.width(9);
		of2.width(9);
		of1<<recgnition<<"  ";*/
		recgnition = allhr[i].cpt_rcgrate(2);
		sumofrate2 = sumofrate2 + recgnition;
		cout<<recgnition<<endl;
		of2<<recgnition<<"   ";
		double rank = 0.0;
		/*rank = allhr[i].get_rank(1);
		sumofrank1 = sumofrank1+rank;
		cout<<"��"<<i<<"����Ȼ�����ʽṹ����Ϊ��"<<rank<<"  ";
		of1.width(9);
		of2.width(9);
		of1<<rank<<"  ";*/
		rank = allhr[i].get_rank(2);
		sumofrank2 = sumofrank2 + rank;
		cout<<rank<<endl;
		of2<<rank<<"   ";
		double relevance = 0.0;
		/*relevance = allhr[i].get_rel(1);
		sumofrlv1 = sumofrlv1 + relevance;
		cout<<"��"<<i<<"�������ʵ������ֵΪ��"<<relevance<<"  ";
		of1.width(9);
		of2.width(9);
		of1<<relevance<<"  ";*/
		//relevance = allhr[i].get_rel(2);
		relevance = allhr[i].get_rel();
		sumofrlv2 = sumofrlv2 + relevance;
		cout<<relevance<<endl;
		of2<<relevance<<"   ";
		double z_score = 0.0;
		/*z_score = allhr[i].get_z(1);
		sumofzs1 = sumofzs1 + z_score;
		cout<<"��"<<i<<"�������ʵ�z_scoreֵΪ��"<<z_score<<"  ";
		of1.width(9);
		of2.width(9);
		of1<<z_score<<"  ";*/
		//z_score = allhr[i].get_z(2);
		z_score = allhr[i].get_z();
		sumofzs2 = sumofzs2 + z_score;
		cout<<z_score<<endl;
		of2<<z_score<<"   ";
		double rmsd = 0.0;
		/*rmsd = allhr[i].count_rmsd(1);
		sumofdrmsd1 = sumofdrmsd1 + rmsd;
		cout<<"��"<<i<<"�������ʵ�drmsd��ֵΪ��"<<rmsd<<"  ";
		of1.width(9);
		of2.width(9);
		of1<<rmsd<<"  ";*/
		rmsd = allhr[i].count_rmsd(2);
		sumofdrmsd2 = sumofdrmsd2 + rmsd;
		cout<<rmsd<<endl;
		of2<<rmsd<<"   ";
		double enrich;
		/*enrich = allhr[i].get_nenrich(1);
		sumofenrich1 = sumofenrich1 + enrich;
		cout<<"��"<<i<<"�������ʵ�denrichment ��ֵΪ��"<<enrich<<"  ";
		of1.width(9);
		of2.width(9);
		of1<<enrich<<endl;*/
		enrich = allhr[i].get_nenrich(2);
		sumofenrich2 = sumofenrich2 + enrich;
		cout<<enrich<<endl;		
		of2<<enrich<<endl;		
	}
	//cout<<"��һ�ַ�������ȷʶ����Ϊ��"<<sumofrate1/allhr.size()<<endl;
	cout<<"�ڶ��ַ�������ȷʶ����Ϊ��"<<sumofrate2/allhr.size()<<endl;
	//cout<<"��һ�ַ�����ƽ����Ϊ��"<<sumofrank1/allhr.size()<<endl;
	cout<<"�ڶ��ַ�����ƽ����Ϊ��"<<sumofrank2/allhr.size()<<endl;
	//cout<<"��һ�ַ�����ƽ�������Ϊ��"<<sumofrlv1/allhr.size()<<endl;
	cout<<"�ڶ��ַ�����ƽ�������Ϊ��"<<sumofrlv2/allhr.size()<<endl;
	//cout<<"��һ�ַ�����ƽ��Z����Ϊ��"<<sumofzs1/allhr.size()<<endl;
	cout<<"�ڶ��ַ�����ƽ��Z����Ϊ��"<<sumofzs2/allhr.size()<<endl;
	//cout<<"��һ�ַ�����ƽ��drmsdΪ��"<<sumofdrmsd1/allhr.size()<<endl;
	cout<<"�ڶ��ַ�����ƽ��drmsdΪ��"<<sumofdrmsd2/allhr.size()<<endl;
	//cout<<"��һ�ַ�����ƽ��Nenrichment��ֵΪ��"<<sumofenrich1/allhr.size()<<endl;
	cout<<"�ڶ��ַ�����ƽ��Nenrichment��ֵΪ��"<<sumofenrich2/allhr.size()<<endl;
	//of1<<"average:"<<endl;
	of2<<"average:"<<endl;

	//of1<<sumofrate1/allhr.size()<<endl;
	of2<<sumofrate2/allhr.size()<<endl;
	//of1<<sumofrank1/allhr.size()<<endl;
	of2<<sumofrank2/allhr.size()<<endl;
	//of1<<sumofrlv1/allhr.size()<<endl;
	of2<<sumofrlv2/allhr.size()<<endl;
	//of1<<sumofzs1/allhr.size()<<endl;
	of2<<sumofzs2/allhr.size()<<endl;
	//of1<<sumofdrmsd1/allhr.size()<<endl;
	of2<<sumofdrmsd2/allhr.size()<<endl;
	//of1<<sumofenrich1/allhr.size()<<endl;
	of2<<sumofenrich2/allhr.size()<<endl;
}