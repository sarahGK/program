#ifndef _RMSD_H
#define _RMSD_H

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include "kabsch.h"

using namespace std;


bool Ev_Rmsd(const string& filename1,const string& filename2, int leibie,double& rmsd);//评估两个蛋白质的RMSD(序列必须一致)

#endif