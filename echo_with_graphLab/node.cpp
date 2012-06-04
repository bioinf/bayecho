#include "node.h"
#include <cmath>
#include <iostream>

map<unsigned int, int> MyNode::number_of_copy; 


MyNode::MyNode(unsigned int read_id1, MMAPReads* readfile1, Options* opt1, int max_seq_len1, double*** loglikelihood1, vector<tr1::tuple<int, int, int> >& hypothesis1){
  max_seq_len = max_seq_len1;
  read_id = read_id1; 
  readfile = readfile1; 
  opt = opt1;
  confMat = new double**[max_seq_len];
  loglikelihood = loglikelihood1;
  for (int i = 0; i < max_seq_len; ++i){
    confMat[i] = new double*[4];
    for (int j = 0; j < 4; ++j){
      confMat[i][j] = new double[4];
    }
  }		
  hypothesis = hypothesis1;
  for (int i = 0; i < 500; ++i){
    histogram.push_back(0);
  }
  if (MyNode::number_of_copy.find(read_id) != MyNode::number_of_copy.end()){
  	cerr<<"FAIL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  }
  nVote = 0;
  MyNode::number_of_copy[read_id] = 1; 
  zero_Mat(max_seq_len, confMat);
}

MyNode::MyNode(const MyNode& node){
  correct_read = node.correct_read;
  qual = node.qual;
  read_id = node.read_id;
  indexs = node.indexs;
  readfile = node.readfile;
  opt = node.opt;
  confMat = node.confMat;
  loglikelihood = node.loglikelihood;
  hypothesis = node.hypothesis;
  max_seq_len = node.max_seq_len;
  histogram = node.histogram;
  hist_readset = node.hist_readset;
  MyNode::number_of_copy[node.read_id] += 1;
  nVote = node.nVote;
}

void zero_Mat(int seq_len, double*** mat) {
  for(int l=0; l<seq_len; l++)
    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
        mat[l][i][j]=0;
}

MyNode::~MyNode(){
  if (MyNode::number_of_copy[read_id] == 1){
    for (int i = 0; i < max_seq_len; ++i){
      for (int j = 0; j < 4; ++j){
        delete[] confMat[i][j];
      }
      delete[] confMat[i];
    }
    delete[] confMat;
  } else {
    MyNode::number_of_copy[read_id] -= 1;
  }
}
			
	
