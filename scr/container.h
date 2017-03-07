class counter{
public:
	vector<vector<vector<vector<int>>>> array; //the count of aisiai-1si-1 appearance
	vector<vector<vector<vector<double>>>> prob;//calculate p(si|aiai-1si-1)
	vector<vector<vector<vector<double>>>> aprob;//calucate P(si|ai-1si-1)*p(si|aiai-1si-1)
	bool flag;				    //true--first definition;false--second definition

	counter(bool f);		      //the initial value for array、prob、aprob 
	void add(int x1,int y1,int x2,int y2);//the array for ai-1、si-1、ai、si + 1
	void comput();                       //calculate prob、aprob according to the value of arrays
	void print(ofstream & of);           //write data into files
};
