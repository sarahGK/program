#ifndef _RMSD_H
#define _RMSD_H

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include "kabsch.h"

using namespace std;


bool Ev_Rmsd(const string& filename1,const string& filename2, int leibie,double& rmsd);//�������������ʵ�RMSD(���б���һ��)

#endif