#define N 0.1
class HR_genel{
private:
	int num_decoys;	//�Ŷ��ṹ����Ŀ
	vector<double> probal1_decoys;//����һ������Ŷ��ṹ�ĸ���ֵ
	vector<double> probal2_decoys;//����һ������Ŷ��ṹ�ĸ���ֵ
	vector<double> rmsd;//�õ����������Ŷ��ṹ��rmsdֵ
	vector<int> prbrank1;//������һ������Ŷ��ṹ����ֵ�����ı�ʾ˳��
	vector<int> prbrank2;//��������������Ŷ��ṹ����ֵ�����ı�ʾ˳��
	double prbvalue1_native;//����һ�������Ȼ�ṹ�ĸ���ֵ
	double prbvalue2_native;//�������������Ȼ�ṹ�ĸ���ֵ
	double recognition_rate;//ʶ����
	int  rank;//��
	double z_score;//z����
	double revalance;//���ϵ��
	double drmsd;
	double n_enrichment;

public:

	string pdbid;//�����ʵ�����
	vector<int> decoys;//alldecoys.txt�м�¼���Ŷ��ṹ�ı�ʾ
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
	double cpt_rcgrate(int no);//no=1��ʾ�Ե�һ����ĸ���巽�����㣬����Ϊ�ڶ��ַ���
	int get_rank(int no);
	double get_z(int no);
	double get_z();//û�в�����ʾ����task3��ֻ�Եڶ��ַ������еļ���
	double get_rel(int no);
	double count_rmsd(int no);
	double get_nenrich(int no);
	int find(int x);
	int find_mark(int x);
	double get_re(int no);
	double get_rel();
};

