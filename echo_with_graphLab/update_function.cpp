#include <iostream>
#include <vector>
#include <tr1/tuple>
#include <fstream>
#include <iostream>
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
#include "update_function.h"
#include "node.h"
#include "edge.h"

gl_types::glshared<double***> CONF_MAT;

void add_to_confMat(graphlab::any& current_confMat, const graphlab::any& added_value){
  int i,j,k;
  pair<vector<int>, double> val = added_value.as<pair< vector<int>, double> >();
  vector<int> dim;
  dim = val.first;
  i = dim[0];
  j = dim[1];
  k = dim[2];
  double value =val.second;
  current_confMat.as<double***>()[i][j][k] += value;
}

void graph_update(gl_types::iscope &scope, gl_types::icallback &scheduler) {
  tr1::shared_ptr<MyNode> vdata= scope.vertex_data();
  if (vdata->get_read_id() >= vdata->opt->read_ed){
    return;
  }
  const string orig_seq = string(vdata->get_readfile()[vdata->get_read_id()]);
  const int seq_len = orig_seq.size();
  string corr_seq(seq_len, 0);
  string qual_seq(seq_len, 0);
  vector<vector<VoteInfo> > Votes(seq_len);
  vector<int> nVotes(seq_len, 0);

  MMAPReads& readfile = vdata->get_readfile();
  unsigned int readid = vdata->get_read_id();
  if(!readfile.isOrig(readid)) {
    for(int i=0; i<seq_len; i++) {
      corr_seq[i] = 'N';
      qual_seq[i] = 33;
    }
    vdata->correct_read = corr_seq;
    vdata->qual = qual_seq;
    vdata->nVote = -1;
    return;
  }
  map<unsigned int, NeighborInfo> readNeighbors;
  map<unsigned int, graphlab::edge_id_t> readid_to_edgeid;
  foreach(graphlab::edge_id_t eid, scope.out_edge_ids()) {
    tr1::shared_ptr<MyNode> neighbor =  scope.neighbor_vertex_data(scope.target(eid));
    readid_to_edgeid.insert(pair<unsigned int, graphlab::edge_id_t>(neighbor->read_id, eid));
    tr1::shared_ptr<Edge> edge= scope.edge_data(eid);
    NeighborInfo nn (edge->get_offset(), edge->get_nerr());
    readNeighbors[neighbor->get_read_id()] = nn;
  }
  NeighborInfo myself(0,0);
  readNeighbors[readid] = myself;


  // Voting.
  // Collect votes.
  for(map<unsigned int, NeighborInfo>::iterator neighbors = readNeighbors.begin(); neighbors!=readNeighbors.end(); neighbors++) {
    const unsigned int& neighborId = neighbors->first;
    const char* neighbor_seq = vdata->get_readfile()[neighborId];
    const int neighbor_seq_len = strlen(neighbor_seq);
    const int overlap = neighbors->second.get_overlap(seq_len, neighbor_seq_len);
    const double log_quality = 0;

    // This check is not redundant, it is used when selecting h/e
    // so that the neighbor computation can be reused.

    if(!neighbors->second.isNeighbor(seq_len, neighbor_seq_len, vdata->opt->K, vdata->opt->h, vdata->opt->e)){
      if (readid_to_edgeid.find(neighbors->first) != readid_to_edgeid.end())
        scope.edge_data(readid_to_edgeid[neighbors->first])->not_neighbor();
      continue;
    }
    const int& st1 = neighbors->second.get_st1();
    const int& st2 = neighbors->second.get_st2();
    const char* overlap_seq = &neighbor_seq[st2];
    const bool orig = readfile.isOrig(neighborId);
    int pos = orig?(st2):(neighbor_seq_len-(st2)-1);
    for(int offset=0; offset<overlap; offset++) {
      int calledb = baseToInt(overlap_seq[offset]);
      if(!isAnyBase(calledb)) {
        if(!orig) calledb = 3-calledb;
        Votes[st1+offset].push_back(VoteInfo(pos, calledb, log_quality, !orig));
        nVotes[st1+offset] += 1;
        // Original read gets to vote twice.
        if(readid==neighborId) {
          int nprior = max(1, (int)floor(sqrt(vdata->opt->cov)));
          if(vdata->opt->h_rate!=0) nprior = 0;
          for(int z=0; z<nprior; z++)
            Votes[st1+offset].push_back(VoteInfo(pos, calledb, log_quality, !orig, true));
        }
      }
      if(orig) pos++; else pos--;
    }
  }


  // Extract most likely sequence and output.
  for(int i=0; i<seq_len; i++) {
    const int orig_base = baseToInt(orig_seq.at(i));
    vector<float> base_loglikelihood(16, 0.0); // Likelihood with prior votes.
    vector<float> base_logquality(16, 0.0);	    // Likelihood without prior votes.

    // perform actual ML estimation

    for(vector<VoteInfo>::iterator pVote=Votes[i].begin(); pVote!=Votes[i].end(); pVote++)
      for(vector<tr1::tuple<int,int,int> >::iterator H=vdata->hypothesis_begin; H!=vdata->hypothesis_end; H++) {
        int& b1 = tr1::get<0>(*H);
        int& b2 = tr1::get<1>(*H);
        int& b1b2 = tr1::get<2>(*H);
        int read_b1 = pVote->reverse_complement?3-b1:b1;
        int read_b2 = pVote->reverse_complement?3-b2:b2;
        if(b1 ==b2) {
          base_loglikelihood[b1b2] += vdata->loglikelihood[pVote->pos][pVote->base][read_b1] + pVote->log_quality;
          if(!pVote->prior)			
            base_logquality[b1b2] += vdata->loglikelihood[pVote->pos][pVote->base][read_b1] + pVote->log_quality;			
        } else  {  
          cerr <<"b1 != b2\n";
        }
      }
    tr1::tuple<int, int, int> mse_base;
    double max_loglikelihood = -numeric_limits<double>::infinity();
    for(vector<tr1::tuple<int,int,int> >::iterator H=vdata->hypothesis_begin; H!=vdata->hypothesis_end; H++) {
      int& b1 = tr1::get<0>(*H);
      int& b2 = tr1::get<1>(*H);
      int& b1b2 = tr1::get<2>(*H);

      if(b1==b2) {
        base_loglikelihood[b1b2] += log(0.25) + log(1.0 - vdata->opt->h_rate);
        base_logquality[b1b2] += log(0.25) + log(1.0 - vdata->opt->h_rate);		    
      } else if(b1<b2) {
        cerr << "b1<b2\n";
      }

      if(base_loglikelihood[b1b2]>max_loglikelihood) {
        max_loglikelihood = base_loglikelihood[b1b2];
        mse_base = *H;
      }
    }

    // Extract ML.
    if(nVotes[i]!=0 && (nVotes[i]<=vdata->opt->max_cov || vdata->opt->max_cov==0) && nVotes[i]>=vdata->opt->min_cov) {
      // accept changes
      corr_seq[i] = intToBase(tr1::get<0>(mse_base), tr1::get<1>(mse_base));	
    } else {
      // Reject correction if received too many/few votes.
      corr_seq[i] = orig_seq.at(i);
      tr1::get<0>(mse_base) = tr1::get<1>(mse_base) = orig_base;
    }

    // Quality score computation.
    double total_prob=0, error_prob=0;
    if(!isAnyBase(tr1::get<0>(mse_base)) && !isAnyBase(tr1::get<0>(mse_base)))
      for(vector<tr1::tuple<int,int,int> >::iterator H=vdata->hypothesis_begin; H!=vdata->hypothesis_end; H++) {
        int& b1 = tr1::get<0>(*H);
        int& b2 = tr1::get<1>(*H);
        int& b1b2 = tr1::get<2>(*H);

        float prob = exp(base_logquality[b1b2]-base_logquality[tr1::get<2>(mse_base)]);
        total_prob += prob;
        if(b1!=tr1::get<0>(mse_base) || b2!=tr1::get<1>(mse_base))
          error_prob += prob;
      }
    error_prob /= total_prob;
    if(error_prob >= 1e-100)
      qual_seq[i] = min(93, max(0, (int)floor(-10.0*log(error_prob) / log(10.0) / max(1.0, (double)nVotes[i])))) + 33;
    else
      qual_seq[i] = 33+93;

    // update confusion matrix
    if(total_prob>0 && !isAnyBase(orig_base))
      for(vector<tr1::tuple<int,int,int> >::iterator H=vdata->hypothesis_begin; H!=vdata->hypothesis_end; H++) {
        int& b1 = tr1::get<0>(*H);
        int& b2 = tr1::get<1>(*H);
        int& b1b2 = tr1::get<2>(*H);

        float prob = exp(base_logquality[b1b2]-base_logquality[tr1::get<2>(mse_base)]);
        vector<int> dim;
        dim.push_back(i);
        dim.push_back(b1);
        dim.push_back(orig_base);
        if (prob > 0){
          CONF_MAT.apply(add_to_confMat,pair<vector<int>, double>(dim, 0.5*prob/total_prob));
          dim[1] = b2;
          CONF_MAT.apply(add_to_confMat, pair<vector<int>, double>(dim, 0.5*prob/total_prob));
        }
      }
  }

  // Update histogram.
  if(vdata->opt->save_stats){
    vdata->nVote = nVotes[0];
  }
  // Output corrected read.
  vdata->correct_read = corr_seq;
  vdata->qual = qual_seq;
}

