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
#include <omp.h>
#include <tbb/concurrent_hash_map.h>
#include <boost/math/special_functions/fpclassify.hpp>

#include "util.hpp"
#include "DNASeq.hpp"
#include "MMAPReads.hpp"
#include "NeighborSet.hpp"

using namespace std;

typedef tbb::concurrent_hash_map<unsigned int, string > conc_map;

unsigned int BlockSize = 10000;

void print_time(const char* str){
  time_t cur_time;
  struct tm * timeinfo;
  time( &cur_time);
  timeinfo = localtime( &cur_time);
  ostringstream str_stream;
  str_stream << asctime(timeinfo) << " " << str << "\n";
  cerr << str_stream.str();
}

struct VoteInfo {
  int pos;
  int base;
  double log_quality;
  bool reverse_complement;
  bool prior;

  VoteInfo(int pos, int base, double log_quality, bool reverse_complement, bool prior=false):pos(pos),base(base),log_quality(log_quality),reverse_complement(reverse_complement),prior(prior) {}
};

void zeroMat(int seq_len, double mat[][4][4]) {
  for(int l=0; l<seq_len; l++)
    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
        mat[l][i][j]=0;
}

void initLoglikelihoodMat(const Options& opt, int max_seq_len, double confMat[][4][4], double loglikelihood[][4][4]) {
  zeroMat(max_seq_len, confMat);
  zeroMat(max_seq_len, loglikelihood);

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
}

void generateHypothesis(bool heterozygous, vector<tr1::tuple<int, int, int> >& hypothesis) {
  hypothesis.clear();

  if(heterozygous) {
    for(int b1 = 0; b1 < 4; b1++)
      for(int b2 = b1; b2 < 4; b2++) {
        hypothesis.push_back(tr1::make_tuple(b1, b2, b1 * 4 + b2));
      }
  } else  {
    for(int b = 0; b < 4; b++)
      hypothesis.push_back(tr1::make_tuple(b, b, b * 4 + b));
  }
}

int main(int argc, char** argv) {
  // Initialize constants.
  Options opt(argc, argv);

  // Initialize reads, votes, and confusion matrix.
  MMAPReads readfile(opt.readFName);
  int max_seq_len = 0;
  for(unsigned int readid=opt.read_st; readid<opt.read_ed; readid++) {
    int seq_len = strlen(readfile[readid]);
    if(seq_len>max_seq_len) max_seq_len = seq_len;
  }

  // candidate hypothesis
  vector<tr1::tuple<int, int, int> > hypothesis;
  generateHypothesis(opt.h_rate>0, hypothesis);

  // Confusion matrix.
  // Error loglikelihood matrix.
  double confMat[max_seq_len][4][4];
  double loglikelihood[max_seq_len][4][4];
  initLoglikelihoodMat(opt, max_seq_len, confMat, loglikelihood);

  // Voting Mechanism.
  // Statitstics.
  vector<int> histogram(500, 0);
  vector<bool> hist_readset_was(opt.read_ed, false);
  // Open output file.
  ostringstream fname, fqualname;
  fname << opt.fpre << "output_" << opt.fsuf << ".txt";
  fqualname << opt.fpre << "quality_" << opt.fsuf << ".txt";
  ofstream fout, fqual;
  fout.open(fname.str().c_str());
  fqual.open(fqualname.str().c_str());

  // initialize NeighborSetLoaders
  NeighborSetLoader  neighborLoader(opt.inputFNames[0]);
  if (opt.inputFNames.size() != 1){
    cerr << "Error: several input files with neighbors\n";
  }
  // Zero out confusion matrix.
  zeroMat(max_seq_len, confMat);
  omp_lock_t conf_mat_lock;
  omp_lock_t hist_lock;
  omp_lock_t out_lock;
  omp_init_lock(&conf_mat_lock);
  omp_init_lock(&hist_lock);
  omp_init_lock(&out_lock);
  unsigned int count_blocks = (opt.read_ed - opt.read_st)/ BlockSize + 1;
  for (unsigned int block_number = 0; block_number < count_blocks; ++ block_number){
    unsigned int begin_block_index = opt.read_st + block_number * BlockSize;
    unsigned int end_block_index = opt.read_st + (block_number+ 1) * BlockSize;
    if (end_block_index > opt.read_ed){
      end_block_index = opt.read_ed;
    }
    cerr << "begin load \n";
    neighborLoader.load(begin_block_index, end_block_index - 1);
    tbb::concurrent_hash_map<unsigned int, string > corr_seqs;
    tbb::concurrent_hash_map<unsigned int, string > corr_quals;

    #pragma omp parallel for  
    for(unsigned int readid = begin_block_index; readid<end_block_index; readid++) {
      const string orig_seq = string(readfile[readid]);
      const int seq_len = orig_seq.size();
      string corr_seq(seq_len, 0);
      string qual_seq(seq_len, 0);
      vector<vector<VoteInfo> > Votes(seq_len);
      vector<int> nVotes(seq_len, 0);

      // If read is reverse complement read, output directly.
      if(!readfile.isOrig(readid)) {
        for(int i=0; i<seq_len; i++) {
          corr_seq[i] = 'N';
          qual_seq[i] = 33;
        }
        conc_map::accessor acc_qual;
        corr_quals.insert(acc_qual, readid);
        acc_qual->second = qual_seq;
        conc_map::accessor acc_seq;
        corr_seqs.insert(acc_seq, readid);
        acc_seq->second = corr_seq;
        continue;
      }
      tr1::shared_ptr<map<unsigned int, NeighborInfo> > newNeighbors = neighborLoader.get(readid);
      // Voting.
      // Collect votes.
      //ostringstream str;
      //str << readid << " "<<newNeighbors->size()<<"\n";
      //cerr << str.str();
      set<unsigned int> my_neighbors;
      for(map<unsigned int, NeighborInfo>::iterator neighbors = newNeighbors->begin();
          neighbors!=newNeighbors->end(); neighbors++) {

        const unsigned int& neighborId = neighbors->first;
        const char* neighbor_seq = readfile[neighborId];
        const int neighbor_seq_len = strlen(neighbor_seq);
        const int overlap = neighbors->second.get_overlap(seq_len, neighbor_seq_len);

        // This check is not redundant, it is used when selecting h/e
        // so that the neighbor computation can be reused.
        if(!neighbors->second.isNeighbor(seq_len, neighbor_seq_len, opt.K, opt.h, opt.e))
          continue;

        if(opt.save_stats)
          my_neighbors.insert(neighborId);

        const int& st1 = neighbors->second.get_st1();
        const int& st2 = neighbors->second.get_st2();
        const char* overlap_seq = &neighbor_seq[st2];
        const bool orig = readfile.isOrig(neighborId);
        int pos = orig?(st2):(neighbor_seq_len-(st2)-1);		
        for(int offset = 0; offset<overlap; offset++) {
          int calledb = baseToInt(overlap_seq[offset]);
          if(!isAnyBase(calledb)) {
            if(!orig) calledb = 3-calledb;
            Votes[st1+offset].push_back(VoteInfo(pos, calledb, 0, !orig));
            nVotes[st1+offset] += 1;
            // Original read gets to vote twice.
            if(readid==neighborId) {
              int nprior = max(1, (int)floor(sqrt(opt.cov)));
              if(opt.h_rate!=0) nprior = 0;
              for(int z=0; z<nprior; z++)
                Votes[st1+offset].push_back(VoteInfo(pos, calledb, 0, !orig, true));
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
          for(vector<tr1::tuple<int,int,int> >::iterator H=hypothesis.begin(); H!=hypothesis.end(); H++) {
            int& b1 = tr1::get<0>(*H);
            int& b2 = tr1::get<1>(*H);
            int& b1b2 = tr1::get<2>(*H);
            int read_b1 = pVote->reverse_complement?3-b1:b1;

            if(b1>b2) {
              cerr<<"error b1> b2\n";
            } else if(b1==b2)  {
              // Homozygous case.
              base_loglikelihood[b1b2] += loglikelihood[pVote->pos][pVote->base][read_b1] + pVote->log_quality;
              if(!pVote->prior)			
                base_logquality[b1b2] += loglikelihood[pVote->pos][pVote->base][read_b1] + pVote->log_quality;			
            } else  { // b1 < b2
              cerr << "error b1 < b2\n";
            }
          }

        tr1::tuple<int, int, int> mse_base;
        double max_loglikelihood = -numeric_limits<double>::infinity();
        for(vector<tr1::tuple<int,int,int> >::iterator H=hypothesis.begin(); H!=hypothesis.end(); H++) {
          int& b1 = tr1::get<0>(*H);
          int& b2 = tr1::get<1>(*H);
          int& b1b2 = tr1::get<2>(*H);

          if(b1==b2) {
            base_loglikelihood[b1b2] += log(0.25) + log(1.0 - opt.h_rate);
            base_logquality[b1b2] += log(0.25) + log(1.0 - opt.h_rate);		    
          } else if(b1<b2) {
            base_loglikelihood[b1b2] += -log(6.0) + log(opt.h_rate);
            base_logquality[b1b2] += -log(6.0) + log(opt.h_rate);		    
          }

          if(base_loglikelihood[b1b2]>max_loglikelihood) {
            max_loglikelihood = base_loglikelihood[b1b2];
            mse_base = *H;
          }
        }
        // Extract ML.
        if(nVotes[i]!=0 && (nVotes[i]<=opt.max_cov || opt.max_cov==0) && nVotes[i]>=opt.min_cov) {
          // accept changes
          corr_seq[i] = intToBase(tr1::get<0>(mse_base), tr1::get<1>(mse_base));		
        } else  {
          // Reject correction if received too many/few votes.
          corr_seq[i] = orig_seq.at(i);
          tr1::get<0>(mse_base) = tr1::get<1>(mse_base) = orig_base;
        }

        // Quality score computation.
        double total_prob=0, error_prob=0;
        if(!isAnyBase(tr1::get<0>(mse_base)) && !isAnyBase(tr1::get<0>(mse_base)))
          for(vector<tr1::tuple<int,int,int> >::iterator H=hypothesis.begin(); H!=hypothesis.end(); H++) {
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
        if( !isinf(total_prob) && !isnan(total_prob) && total_prob>0 && !isAnyBase(orig_base))
          for(vector<tr1::tuple<int,int,int> >::iterator H=hypothesis.begin(); H!=hypothesis.end(); H++) {
            int& b1 = tr1::get<0>(*H);
            int& b2 = tr1::get<1>(*H);
            int& b1b2 = tr1::get<2>(*H);

            float prob = exp(base_logquality[b1b2]-base_logquality[tr1::get<2>(mse_base)]);
            if(readfile.isOrig(readid)) {
              if (prob > 0){
                omp_set_lock(&conf_mat_lock);
                confMat[i][b1][orig_base] += 0.5*prob/total_prob;
                confMat[i][b2][orig_base] += 0.5*prob/total_prob;
                omp_unset_lock(&conf_mat_lock);
              }
            } else  {
              cerr << "not origin\n";
            }
          }
      }
      omp_set_lock(&hist_lock);
      // Update histogram.
      if(opt.save_stats){
        size_t max_vote = *max_element(nVotes.begin(), nVotes.end());
        if(max_vote>histogram.size()){
          histogram.resize(max_vote, 0);
        }
        if(!hist_readset_was[readid]){
          // All neighbors will not vote in the future and only beginning of the read gets to vote to improve independence.
          for (set<unsigned int>::iterator it = my_neighbors.begin(); it != my_neighbors.end(); ++it){
            if (*it < opt.read_ed)
              hist_readset_was[*it] = true;
          }
          if(readfile.isOrig(readid))
            histogram[nVotes[0]]+=1;
          else
            histogram[nVotes[seq_len-1]]+=1;
        }
      }
      omp_unset_lock(&hist_lock);
      // Output corrected read.
      conc_map::accessor acc_qual;
      corr_quals.insert(acc_qual, readid);
      acc_qual->second = qual_seq;
      conc_map::accessor acc_seq;
      corr_seqs.insert(acc_seq, readid);
      acc_seq->second = corr_seq;
      newNeighbors->clear();
      newNeighbors.reset();
    }
    for(unsigned int readid = begin_block_index; readid < end_block_index; readid++){
      conc_map::const_accessor ac_qual;
      corr_quals.find(ac_qual, readid);
      conc_map::const_accessor ac_seq;
      corr_seqs.find(ac_seq, readid);
      fout << ac_seq->second << endl;
      fqual << ac_qual->second << endl;
    }
    corr_quals.clear();
    corr_seqs.clear();
  }

  fout.close();
  fqual.close();
  // Output other stats.
  // Output histogram.
  if(opt.save_stats) {
    fname.str(string());
    fname << opt.fpre << "histogram_" << opt.fsuf << ".txt";
    fout.open(fname.str().c_str());
    for(size_t i=0; i<histogram.size(); i++)
      fout << histogram[i] << ' ';
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
  return 0;
}
