#ifndef _RD_FASTA_H
#define _RD_FASTA_H
#define MAX_LINE 8000
#include <fstream>
#include <string>
#include <vector>
using namespace std;

bool ReadFastaSeq(ifstream& file,string &seq);
bool ReadFastaSeq(ifstream& file,string &seq,string& miaoshu);
bool ReadFastaSeq(ifstream& file,string &seq,char& cc);
bool ReadFastaSeq(ifstream& file,string &seq,string& miaoshu,char& cc);
bool WriteFastaSeq(string& file,string &seq);
bool WriteFastaSeq(string& filename,string &seq,string& miaoshu);
bool WriteFastaSeq(ofstream& file,string &seql);
bool WriteFastaSeq(ofstream& file,string &seq,string& miaoshu);
#endif