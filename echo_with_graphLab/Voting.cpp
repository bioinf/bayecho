	#include <iostream>
#include <vector>
#include <tr1/tuple>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <utility>
#include <algorithm>
#include <tr1/memory>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>

#include "util.hpp"
#include "DNASeq.hpp"
#include "MMAPReads.hpp"
#include "NeighborSet.hpp"

#include "node.h"
#include "edge.h"
#include "update_function.h"
#include"result_histogram.h"

using namespace std;


typedef std::map<unsigned int, std::tr1::shared_ptr<std::map<unsigned int, NeighborInfo> > >  NeighborMap;

gl_types::glshared<double> RED_PROPORTION;

void deleteMatrix(double*** matrix, int len){
  for (int i = 0; i < len; ++i){
    for (int j = 0; j < 4; ++j){
      delete[] matrix[i][j];
    }
    delete[] matrix[i];
  }
  delete[] matrix;  
}

double*** initLoglikelihood_Mat(const Options& opt, int max_seq_len) {
  double*** confMat = new double**[max_seq_len];
  double*** loglikelihood = new double**[max_seq_len];
  for (int i = 0; i < max_seq_len; ++i){
    confMat[i] = new double*[4];
    loglikelihood[i] = new double*[4];
    for (int j = 0; j < 4; ++j){
      confMat[i][j] = new double[4];
      loglikelihood[i][j] = new double[4];
    }
  }		
  zero_Mat(max_seq_len, confMat);
  zero_Mat(max_seq_len, loglikelihood);

  // Initialize error loglikelihood matrix.
  if(opt.confMatFName==NULL) {
    // If no confusion matrix supplied, use majority weight instead.
    for(int pos=0; pos<max_seq_len; pos++)
      for(int b1=0; b1<4; b1++)
        for(int b2=0; b2<4; b2++)
          if(b1==b2)
            loglikelihood[pos][b1][b2]=log(0.99);
          else
            loglikelihood[pos][b1][b2]=log(0.01/3);
  } else  {
    // Otherwise, load confusion matrix.
    ifstream fin(opt.confMatFName);
    vector<vector<int> > tot(max_seq_len);
    int file_len;
    fin >> file_len;
    for(int pos=0; pos<file_len; pos++) {
      tot[pos].resize(4);
      for(int trueb=0; trueb<4; trueb++)
        for(int calledb=0; calledb<4; calledb++) {
          fin >> confMat[pos][trueb][calledb];
          confMat[pos][trueb][calledb]+=1.0;		    
          tot[pos][calledb] += confMat[pos][trueb][calledb];
        }
      }
      fin.close();

      // Normalize and take the log. 
      for(int pos=0; pos<file_len; pos++)
        for(int trueb=0; trueb<4; trueb++)
          for(int calledb=0; calledb<4; calledb++)
            loglikelihood[pos][calledb][trueb] = log(confMat[pos][trueb][calledb]) - log(tot[pos][calledb]);

      // Fill the gap. 
      for(int pos=file_len; pos<max_seq_len; pos++)
        for(int b1=0; b1<4; b1++)
          for(int b2=0; b2<4; b2++)
            loglikelihood[pos][b1][b2] = loglikelihood[file_len-1][b1][b2];
  }  
  deleteMatrix(confMat, max_seq_len);
  return loglikelihood;  
}

void generateHypothesis(bool heterozygous, vector<tr1::tuple<int, int, int> >& hypothesis) {
  hypothesis.clear();
  if(heterozygous) {
    for(int b1=0; b1<4; b1++)
      for(int b2=b1; b2<4; b2++) {
         hypothesis.push_back(tr1::make_tuple(b1, b2, b1*4+b2));
      }
  } else  {
    for(int b=0; b<4; b++)
      hypothesis.push_back(tr1::make_tuple(b, b, b*4+b));
  }
}

void saveStat(const Options& opt, vector<int>& histogram, int max_seq_len, double*** confMat){
  if(opt.save_stats) {
    // Output histogram.
    ostringstream fname;
    fname << opt.fpre << "histogram_" << opt.fsuf << ".txt";
    ofstream fout;
    fout.open(fname.str().c_str());
    for(size_t i=0; i<histogram.size(); i++){
      fout << histogram[i] << ' ';
    }
    fout.close();
    
    // Output confusion matrix.
    fname.str(string());
    fname << opt.fpre << "confmat_" << opt.fsuf << ".txt";
    fout.open(fname.str().c_str());
    fout << max_seq_len << endl;
    for(int pos=0; pos<max_seq_len; pos++) {
      for(int b1=0; b1<4; b1++) {
        for(int b2=0; b2<4; b2++)
          fout << confMat[pos][b1][b2] << ' ';
        fout << endl;
      }
      fout << endl;
    }
    fout.close();
  }
}


int main(int argc, char** argv) {
  // Initialize constants.
  Options opt(argc, argv);
  // Initialize reads, votes, and confusion matrix.
  MMAPReads* readfile1 = new MMAPReads(opt.readFName);
  MMAPReads readfile = *readfile1;
    
  int max_seq_len = 0;
  int min_read_id = 10000000000;
  int max_read_id = 0;
  for(unsigned int readid=opt.read_st; readid<opt.read_ed; readid++) {
    int seq_len = strlen(readfile[readid]);
    if(seq_len>max_seq_len) 
      max_seq_len = seq_len;
    if (readid < min_read_id){
      min_read_id = readid;
    } 
    if (readid > max_read_id){
      max_read_id = readid;
    }
  }
  double*** loglikelihood = initLoglikelihood_Mat(opt, max_seq_len);
  // candidate hypothesis
  vector<tr1::tuple<int, int, int> > hypothesis;
  generateHypothesis(opt.h_rate > 0, hypothesis);

  // Voting Mechanism.
    
  // Open output file.
  ostringstream fname, fqualname;
  fname << opt.fpre << "output_" << opt.fsuf << ".txt";
  fqualname << opt.fpre << "quality_" << opt.fsuf << ".txt";
  ofstream fout, fqual;
  fout.open(fname.str().c_str());
  fqual.open(fqualname.str().c_str());
 
  // initialize NeighborSetLoaders
  vector<tr1::shared_ptr<NeighborSetLoader> > neighborLoader;
  for(size_t fiter=0; fiter<opt.inputFNames.size(); fiter++)
    neighborLoader.push_back(tr1::shared_ptr<NeighborSetLoader>(new NeighborSetLoader(opt.inputFNames[fiter])));
  //create graph
  gl_types::core core;
  map<int, graphlab::vertex_id_t> readid_to_vertexid; 
  graph_type&  graph= core.graph(); 
  
  for (unsigned int readid=opt.read_st; readid<opt.read_ed; readid++){
    graphlab::vertex_id_t  id = graph.add_vertex(MyNode(readid, readfile1, &opt, max_seq_len, loglikelihood, hypothesis));
    readid_to_vertexid.insert(pair<int, graphlab::vertex_id_t>(readid, id));
  }
  for (int neighb_i = 0; neighb_i < neighborLoader.size(); neighb_i++){
    NeighborMap neighbors = neighborLoader[neighb_i]->get_map(min_read_id, max_read_id);	    
    for ( NeighborMap::iterator it = neighbors.begin(); it != neighbors.end(); ++it){
      if (readid_to_vertexid.find(it->first) == readid_to_vertexid.end()){
        graphlab::vertex_id_t  id = graph.add_vertex(MyNode(it->first, readfile1, &opt, max_seq_len, loglikelihood, hypothesis));
        readid_to_vertexid.insert(pair<int, graphlab::vertex_id_t>(it->first, id));
      }
    }

    set<pair<unsigned int, unsigned int> > edges;
    for ( NeighborMap::iterator it = neighbors.begin(); it != neighbors.end(); ++it){
     for(map<unsigned int, NeighborInfo>::iterator it1 = it->second->begin(); it1 != it->second->end(); ++it1){	
        if (it->first != it1->first){
        	if (readid_to_vertexid[it->first] == readid_to_vertexid[it1->first]){
        	} else {
            if (edges.find(pair<unsigned int, unsigned int>(readid_to_vertexid[it->first], readid_to_vertexid[it1->first])) == edges.end()){
          		graph.add_edge(readid_to_vertexid[it->first], readid_to_vertexid[it1->first], *(new Edge(it1->second.get_offset(), it1->second.get_nerr())));
              edges.insert(pair<unsigned int, unsigned int>(readid_to_vertexid[it->first], readid_to_vertexid[it1->first]));
            }
        	}
	}
      }
    }
  }	
  graph.finalize(); 
  for (graphlab::vertex_id_t vid = 0; vid < graph.num_vertices(); ++vid) {   
    core.add_task(vid, graph_update, 100.);
  }
  core.start();

  //reduce function  
  tr1::shared_ptr<MyResult> result(new MyResult(max_seq_len, opt));
 
  for(unsigned int readid=opt.read_st; readid<opt.read_ed; readid++) {
  //for (graphlab::vertex_id_t vid = 0; vid < graph.num_vertices(); ++vid) { 
    graphlab::vertex_id_t vid = readid_to_vertexid[readid]; 
    //cout << "vid "<< vid << " readid " << readid;
    fout << graph.vertex_data(vid).correct_read<< endl;
    fqual << graph.vertex_data(vid).qual << endl;
    
    if ((graph.vertex_data(vid).correct_read.size() != graph.vertex_data(vid).qual.size()) || graph.vertex_data(vid).qual.size()!=  strlen(readfile[readid]) ){
      cerr <<"Error "<<readid<< " cor_size "<<graph.vertex_data(vid).correct_read.size() << " qual size " <<graph.vertex_data(vid).qual.size() << " str size \n" <<readfile[readid]<<"\n" ;
      cerr <<graph.vertex_data(vid).read_id<<" vid"<< vid <<"\n";
    }


    if(!readfile.isOrig(graph.vertex_data(vid).read_id))
      continue;
    for (int i1 = 0; i1 < max_seq_len; ++i1){
      for (int i2 = 0 ; i2 < 4; ++i2){
        for (int i3 = 0; i3 < 4; ++i3){
          result->confMat[i1][i2][i3] += graph.vertex_data(vid).confMat[i1][i2][i3];
        }
      }
    }

    if(result-> hist_read_ids.find(graph.vertex_data(vid).read_id)==result-> hist_read_ids.end()){
      for (set<unsigned int>::iterator aa = graph.vertex_data(vid).my_neighbors.begin(); aa !=graph.vertex_data(vid).my_neighbors.end(); ++aa){
        result->hist_read_ids.insert(*aa);
      }
      result->histogram[graph.vertex_data(vid).nVote]+=1;	
    }
  }
  fout.close();
  fqual.close(); 
  saveStat(opt, result->histogram, max_seq_len, result->confMat);
  deleteMatrix(loglikelihood, max_seq_len);
  return 0;
}
