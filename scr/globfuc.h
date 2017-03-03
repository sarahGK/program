const char rname_one_rank20[21]="ACDEFGHIKLMNPQRSTVWY";

int BinarySearchChar(char c)
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
}