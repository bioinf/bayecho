#ifndef _MY_NODE_H
#define _MY_NODE_H

#include <string>
#include <set>
#include <vector>
#include "util.hpp"
#include <tr1/tuple>
#include "DNASeq.hpp"
#include "MMAPReads.hpp"
#include "NeighborSet.hpp"
#include <map>
using namespace std;

void zero_Mat(int seq_len, double*** mat);




class MyNode{
	public:
		string correct_read;
		string qual;
		unsigned int read_id;
		set<int> indexs;
		MMAPReads* readfile;
		Options* opt;
		double*** confMat;
    		double*** loglikelihood;
    		vector<tr1::tuple<int, int, int> > hypothesis;
    		int max_seq_len;
    		vector<int> histogram;
    		set<unsigned int> hist_readset;
    		static map<unsigned int, int> number_of_copy;
    		int nVote;
    		set<unsigned int> my_neighbors;
    		map<unsigned int, char> correct_read_chars; 
    		
	public:		
		MyNode(unsigned int read_id1, MMAPReads* readfile1, Options* opt1, int max_seq_len1, double*** loglikelihood1, vector<tr1::tuple<int, int, int> >& hypothesis);
		MyNode(const MyNode& node);
		unsigned int get_read_id(){return read_id;};
		MMAPReads& get_readfile(){ return *readfile;};
		~MyNode(); 
		
};


#endif
