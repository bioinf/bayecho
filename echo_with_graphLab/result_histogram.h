#ifndef _RESULT_HISTOGRAM_
#define _RESULT_HISTOGRAM_
#include<vector>
#include<string>
#include<iostream>
#include<set>
using namespace std;
struct MyResult{
  double*** confMat;
  int max_seq_len;
  vector<int> histogram;
  set<string> correct_reads;
  set<string> quals;
  set<unsigned int> hist_read_ids;

  MyResult(int max_seq_len1, const Options& opt1){
    max_seq_len = max_seq_len1;
    confMat = new double**[max_seq_len1];
    for (int i = 0; i < max_seq_len1; ++i){
      confMat[i] = new double*[4];
      for (int j = 0; j < 4; ++j){
        confMat[i][j] = new double[4];
        for (int k = 0; k < 4; ++k){
          confMat[i][j][k] = 0;
        }
      }
    }	
    for (int i = 0; i < 500; ++i){
      histogram.push_back(0);
    }
  }

  MyResult(){
  }

  MyResult(const MyResult& result){
    confMat = result.confMat;
    max_seq_len = result.max_seq_len;
    histogram = result.histogram;
    correct_reads = result.correct_reads;
    quals = result.quals;
  }	


  ~MyResult(){
    for (int i = 0; i < max_seq_len; ++i){
      for (int j = 0; j < 4; ++j){
        delete[] confMat[i][j];
      }
      delete[] confMat[i];
    }
    delete[] confMat;
  }
};
#endif
