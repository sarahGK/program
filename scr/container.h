class counter{
public:
	vector<vector<vector<vector<int>>>> array;//ͳ��aisiai-1si-1���ֵĴ���
	vector<vector<vector<vector<double>>>> prob;//��¼p(si|aiai-1si-1)
	vector<vector<vector<vector<double>>>> aprob;//��¼P(si|ai-1si-1)*p(si|aiai-1si-1)
	bool flag;//������ĸ���巽��һ(true)��(false)

	counter(bool f);//��array�ĳ�ʼֵ��Ϊ1����һarray��prob��aprob��ά��
	void add(int x1,int y1,int x2,int y2);//�����ֱ��Ӧai-1��si-1��ai��si�������Ӧ��array��һ
	void comput();//��array��ֵ������prob��aprob��ֵ
	void print(ofstream & of);//�����ָ���ļ�
};