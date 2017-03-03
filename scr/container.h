class counter{
public:
	vector<vector<vector<vector<int>>>> array;//统计aisiai-1si-1出现的次数
	vector<vector<vector<vector<double>>>> prob;//记录p(si|aiai-1si-1)
	vector<vector<vector<vector<double>>>> aprob;//记录P(si|ai-1si-1)*p(si|aiai-1si-1)
	bool flag;//区别字母表定义方案一(true)二(false)

	counter(bool f);//设array的初始值都为1，第一array、prob、aprob的维数
	void add(int x1,int y1,int x2,int y2);//参数分别对应ai-1、si-1、ai、si，将其对应的array加一
	void comput();//由array的值来计算prob、aprob的值
	void print(ofstream & of);//输出到指定文件
};