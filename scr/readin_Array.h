//const char rname_one_rank20[21] = "ACDEFGHIKLMNPQRSTVWY";
double proposal1[20][1296];//task2中第一种方案定义字母表时的概率值
double proposal2[20][15];//task2中第一种方案定义字母表时的概率值

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



void initial(){//读入task2的概率数组
	ifstream datain;
	datain.open("E:\\talk\\biology\\task2\\probalistics.txt",ios::in);
	//datain.open("G:\\task3\\result\\probability2.txt",ios::in);
	if(!datain.is_open()){
		cout<<"can't open the file!"<<endl;
		return;
	}
	cout<<"正在读取文件："<<"E:\\talk\\biology\\task2\\probalistics.txt"<<endl;
	//cout<<"正在读取文件："<<"G:\\task3\\result\\probability2.txt"<<endl;
	for(int i=0;i<20;i++)
		for(int j=0;j<1296;j++)
			datain>>proposal1[i][j];
	for(int i=0;i<20;i++)
		for(int j =0;j<15;j++)
			datain>>proposal2[i][j];
}
/*
*以读入的方式打开指定的文件
*/
bool openfile(string filename,ifstream & datain){
	
	datain.open (filename.c_str(),ios::in);
	if(!datain.is_open()){
		printf("无法打开文件：%s\n",filename.c_str());
		return false;
	}
	else{
		cout<<"正在读取文件："<<filename.c_str()<<endl;
		return true;
	}
}
/*
*task2中对方案一每个字母出现的概率值取对数
*/
double log_probal1(char ch,int x){
	
	int i = BinarySearchChar(ch);
	if(i!=-1&&x>=0&&x<1296)	return log(proposal1[i][x]);
	else return 0;
}
/*
*task2中对方案二每个字母出现的概率值取对数
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
*task3中对方案二每个字母出现的概率值取对数
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
*将hr中文件的字母表信息读入allhr中
*/
void prod_hrall(vector<HR_genel> & allhr,counter & c){
    
	ifstream name;
	string nameline;
	int i = 1;
	
	if(openfile("E:\\talk\\biology\\task1\\alldecoys.txt",name)){//按照记录扰动结构文件名的顺序
			HR_genel hr;
			bool flag1=false;//是否读入天然结构
			bool flag2=false;//是否读入500个扰动结构
			bool flag3=false;//是否读入扰动结构的标示
		while(!name.eof()){
		
			if(i==1){//跳过文件第一行
				i++;
				getline(name,nameline);			
				continue;
			}
			if(i%2==0){//确定读取两面角的文件名
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
				
					if(dihedralline.length()==17){//当前读取的是天然结构的信息
						flag1 = true;
						cout<<"正在计算天然结构的概率值．．．"<<endl;
						getline(dihedraldata,dihedralline);
						double sum1,sum2;
						sum1= sum2 = 0.0;
						bool first = true;//标示是否是第一个氨基酸分子
						int x;//前氨基酸的两面角字母
						char c1;//记录前氨基酸的字母表示
						while(dihedralline.length()>20){
						char ch = dihedralline.substr(8,1).c_str()[0];
						//int x = atoi(dihedralline.substr(30,10).c_str());
						int y = atoi(dihedralline.substr(40,12).c_str());
						//sum1 = sum1 + log_probal1(ch,x);
						if(first){
							first = false;
							sum2 = sum2 + log_probal2(ch,y);//第一个氨基酸分子，按照task2计算概率值
						}else sum2 = sum2 + log_probal2(c1,x,ch,y,c);//其余氨基酸按照task3计算概率值
						c1 = ch;//当前值为下一个计算的前氨基酸的字母表示
						x = y;//当前值为下一个计算的前氨基酸两面角的字母表示
						getline(dihedraldata,dihedralline);
						}
						//hr.set_prbvalue(sum1,sum2);
						hr.set_prbvalue(sum2);
			
					}
					
					if(flag1){//读取了天然结构的信息后，读取扰动结构的信息
					
						cout<<"正在计算扰动结构的概率值．．．"<<endl;
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
		
			else{//alldecoys.txt的奇数行记录了取扰动结构的标示
				cout<<"正在读取扰动结构的标示．．．"<<endl;
				int m[500] ;
				for(int k=0;k<500;k++){	
				name>>m[k];
				hr.set_mark(m[k]);
				}
				getline(name,nameline);
				i++;
				flag3 = true;
			}
			
			if(flag1&&flag2&&flag3)	{//当每一个蛋白质的天然结构、扰动结构以及标示都读取完毕，将当前的hr添加到数组中
				allhr.push_back(hr);	
				hr.~HR_genel();
			}
			cout<<allhr.size()<<endl;
			if(!flag3||!flag1||!flag2) continue;
			//HR_genel hr;
			flag1=flag2=flag3=false;//处理下一个蛋白质分子
		}
	}
}

/*
*读入扰动结构的rmsd值
*/
void rmsd_allhr(vector<HR_genel> & allhr){
	for(int i=0;i<allhr.size();i++){
		allhr[i].set_prbrank();//对每一个hr的概率值排序
		 ifstream rmsdfile;
		 string filename="G:\\hr\\hrdecoys_selected\\";
		 filename.append(allhr[i].pdbid.c_str(),5);
		 filename.append("\\",1);
		 filename.append(allhr[i].pdbid.c_str(),5);
		 filename.append(".rmsd_recal",11);
		 if(openfile(filename,rmsdfile)){
				cout<<"正在读入扰动结构的rmsd值．．．"<<endl;
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
*计算6个评价指标值及所有蛋白质的平均值,并输出到指定文件
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
	of1<<"pdbid"<<"   "<<"识别率"<<"   "<<"秩"<<"   "<<"相关性"<<"   "<<"z分数"<<"   "<<"drmsd"<<"   "<<"Nenrichment"<<endl;
	of2<<"pdbid"<<"   "<<"识别率"<<"   "<<"秩"<<"   "<<"相关性"<<"   "<<"z分数"<<"   "<<"drmsd"<<"   "<<"Nenrichment"<<endl;
	of1.width(9);
	of2.width(9);
	for(int i=0;i<allhr.size();i++){
		of1<<allhr[i].pdbid<<"   ";
		of2<<allhr[i].pdbid<<"   ";
		double recgnition = 0.0;
		/*recgnition = allhr[i].cpt_rcgrate(1);
		sumofrate1 = sumofrate1 + recgnition;
		cout<<"第"<<i<<"个蛋白质结构的识别率为："<<recgnition<<"  ";
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
		cout<<"第"<<i<<"个天然蛋白质结构的秩为："<<rank<<"  ";
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
		cout<<"第"<<i<<"个蛋白质的相关性值为："<<relevance<<"  ";
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
		cout<<"第"<<i<<"个蛋白质的z_score值为："<<z_score<<"  ";
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
		cout<<"第"<<i<<"个蛋白质的drmsd的值为："<<rmsd<<"  ";
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
		cout<<"第"<<i<<"个蛋白质的denrichment 的值为："<<enrich<<"  ";
		of1.width(9);
		of2.width(9);
		of1<<enrich<<endl;*/
		enrich = allhr[i].get_nenrich(2);
		sumofenrich2 = sumofenrich2 + enrich;
		cout<<enrich<<endl;		
		of2<<enrich<<endl;		
	}
	//cout<<"第一种方案的正确识别率为："<<sumofrate1/allhr.size()<<endl;
	cout<<"第二种方案的正确识别率为："<<sumofrate2/allhr.size()<<endl;
	//cout<<"第一种方案的平均秩为："<<sumofrank1/allhr.size()<<endl;
	cout<<"第二种方案的平均秩为："<<sumofrank2/allhr.size()<<endl;
	//cout<<"第一种方案的平均相关性为："<<sumofrlv1/allhr.size()<<endl;
	cout<<"第二种方案的平均相关性为："<<sumofrlv2/allhr.size()<<endl;
	//cout<<"第一种方案的平均Z分数为："<<sumofzs1/allhr.size()<<endl;
	cout<<"第二种方案的平均Z分数为："<<sumofzs2/allhr.size()<<endl;
	//cout<<"第一种方案的平均drmsd为："<<sumofdrmsd1/allhr.size()<<endl;
	cout<<"第二种方案的平均drmsd为："<<sumofdrmsd2/allhr.size()<<endl;
	//cout<<"第一种方案的平均Nenrichment的值为："<<sumofenrich1/allhr.size()<<endl;
	cout<<"第二种方案的平均Nenrichment的值为："<<sumofenrich2/allhr.size()<<endl;
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