#include "node.h"
#include <cmath>
#include <iostream>

MyNode::MyNode(unsigned int read_id1, MMAPReads* readfile1, Options* opt1, int max_seq_len1, double*** loglikelihood1, vector<tr1::tuple<int, int, int> >& hypothesis1){
  max_seq_len = max_seq_len1;
  read_id = read_id1; 
  readfile = readfile1; 
  opt = opt1;
  loglikelihood = loglikelihood1;
  hypothesis_begin = hypothesis1.begin();
  hypothesis_end = hypothesis1.end();
  nVote = 0;
}

void zero_Mat(int seq_len, double*** mat) {
  for(int l=0; l<seq_len; l++)
    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
        mat[l][i][j]=0;
}

MyNode::~MyNode(){
}


