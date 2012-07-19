#ifndef __OPTION_H__
#define __OPTION_H__

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

struct Options {
  // DNA read file.
  char *readFName;

  // Other input file names.
  std::vector<char*> inputFNames;

  // Confusion matrix file name.
  char *confMatFName;

  // Clustering parameters.
  int K;
  int h;
  float e;
  unsigned int nhash;
  unsigned int ihash_st;
  unsigned int ihash_ed;    
  int blocksize;

  // Coverage upper bound.
  int max_cov;
  int min_cov;
  int cov;

  // Chunk of reads [read_st, read_ed).
  unsigned int read_st;
  unsigned int read_ed;
  unsigned int read_st2;
  unsigned int read_ed2;

  // File suffix.
  char *fpre;    
  char *fsuf;

  // Output statistics?
  bool save_stats;

  // Heterozygous rate.
  float h_rate;

  Options(int argc, char** argv): readFName(NULL),confMatFName(NULL),K(10),h(20),e(0.3),nhash(1),ihash_st(0),ihash_ed(1),blocksize(1000),max_cov(0),min_cov(0),cov(0),read_st(0),read_ed(0),read_st2(0),read_ed2(0),fpre(NULL),fsuf(NULL),save_stats(false),h_rate(0.0) {
    for(int i=1; i<argc-1; i+=2) {
      if(strcmp(argv[i], "-K")==0)
        K = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-h")==0)
        h = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-e")==0)
        e = atof(argv[i+1]);
      else if(strcmp(argv[i], "-max_cov")==0)
        max_cov = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-min_cov")==0)
        min_cov = atoi(argv[i+1]);	    
      else if(strcmp(argv[i], "-cov")==0)
        cov = atoi(argv[i+1]);	    
      else if(strcmp(argv[i], "-nhash")==0)
        nhash = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-ihash_st")==0)
        ihash_st = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-ihash_ed")==0)
        ihash_ed = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-blocksize")==0)
        blocksize = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-st")==0)
        read_st = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-ed")==0)
        read_ed = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-st2")==0)
        read_st2 = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-ed2")==0)
        read_ed2 = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-fpre")==0)
        fpre = argv[i+1];
      else if(strcmp(argv[i], "-fsuf")==0)
        fsuf = argv[i+1];
      else if(strcmp(argv[i], "-read")==0)
        readFName = argv[i+1];
      else if(strcmp(argv[i], "-input")==0)
        inputFNames.push_back(argv[i+1]);
      else if(strcmp(argv[i], "-confMat")==0)
        confMatFName = argv[i+1];
      else if(strcmp(argv[i], "-save_stats")==0)
        save_stats = atoi(argv[i+1]);
      else if(strcmp(argv[i], "-h_rate")==0)
        h_rate = atof(argv[i+1]);
      else {
        std::cerr << "Unrecognized option: " << argv[i] << std::endl;
        exit(1);
      }
    }
  }
};

#endif
