#define N 0.1
class HR_genel{
private:
	int num_decoys;	//扰动结构的数目
	vector<double> probal1_decoys;//方案一计算的扰动结构的概率值
	vector<double> probal2_decoys;//方案一计算的扰动结构的概率值
	vector<double> rmsd;//该蛋白质所有扰动结构的rmsd值
	vector<int> prbrank1;//按方案一计算的扰动结构概率值排序后的标示顺序
	vector<int> prbrank2;//按方案二计算的扰动结构概率值排序后的标示顺序
	double prbvalue1_native;//方案一计算的天然结构的概率值
	double prbvalue2_native;//方案二计算的天然结构的概率值
	double recognition_rate;//识别率
	int  rank;//秩
	double z_score;//z分数
	double revalance;//相关系数
	double drmsd;
	double n_enrichment;

public:

	string pdbid;//蛋白质的名称
	vector<int> decoys;//alldecoys.txt中记录的扰动结构的标示
	HR_genel();
	~HR_genel();
	void set_ID(string st){pdbid = st;}
	void set_mark(int x){this->decoys.push_back(x);}
	void set_pdecoys(double d1,double d2){probal1_decoys.push_back(d1);probal2_decoys.push_back(d2);}
	void set_pdecoys(double d){probal2_decoys.push_back(d);}
	void set_rmsd(double d){rmsd.push_back(d);}
	void set_prbvalue(double d1,double d2){prbvalue1_native = d1;prbvalue2_native = d2;}
	void set_prbvalue(double d){prbvalue2_native = d;}
	void set_prbrank();
	double cpt_rcgrate(int no);//no=1表示对第一种字母定义方案计算，否则为第二种方案
	int get_rank(int no);
	double get_z(int no);
	double get_z();//没有参数表示，在task3中只对第二种方案进行的计算
	double get_rel(int no);
	double count_rmsd(int no);
	double get_nenrich(int no);
	int find(int x);
	int find_mark(int x);
	double get_re(int no);
	double get_rel();
};

