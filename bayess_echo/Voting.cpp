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
#include <time.h>


#include "util.hpp"
#include "DNASeq.hpp"
#include "MMAPReads.hpp"
#include "NeighborSet.hpp"
#include "brother.h"
#include "hammer.h"

using namespace std;

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

void initLoglikelihoodMat(const Options& opt, int max_seq_len, double confMat[][4][4], double loglikelihood[][4][4], double*** confMat1) {
  zeroMat(max_seq_len, confMat);
  zeroMat(max_seq_len, loglikelihood);

  // Initialize error loglikelihood matrix.
  if(opt.confMatFName==NULL) {
    // If no confusion matrix supplied, use majority weight instead.
    for(int pos=0; pos<max_seq_len; pos++)
      for(int b1=0; b1<4; b1++)
        for(int b2=0; b2<4; b2++)
          if(b1==b2){
            confMat1[pos][b1][b2] = 0.99;
            loglikelihood[pos][b1][b2]=log(0.99);
          } else{
            confMat1[pos][b1][b2] = 0.01/3;
            loglikelihood[pos][b1][b2]=log(0.01/3);
          }
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
        for(int calledb=0; calledb<4; calledb++){
          loglikelihood[pos][calledb][trueb] = log(confMat[pos][trueb][calledb]) - log(tot[pos][calledb]);
          confMat1[pos][trueb][calledb] = confMat[pos][trueb][calledb]/tot[pos][calledb];
        }

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
    for(int b1=0; b1<4; b1++)
      for(int b2=b1; b2<4; b2++) {
        hypothesis.push_back(tr1::make_tuple(b1, b2, b1*4+b2));
      }
  } else  {
    for(int b=0; b<4; b++)
      hypothesis.push_back(tr1::make_tuple(b, b, b*4+b));
  }
}
void print_time(const char* str){
  time_t cur_time;
  struct tm * timeinfo;
  time( &cur_time);
  timeinfo = localtime( &cur_time); 
  cerr << str << " "  << asctime(timeinfo) << "\n";
}

int main(int argc, char** argv) {
  print_time("Begin voting");
  // Initialize constants.
  Options opt(argc, argv);
  
  // Initialize reads, votes, and confusion matrix.
  MMAPReads readfile(opt.readFName);
  print_time("constructed read file");
  int max_seq_len = 0;
  for(unsigned int readid=opt.read_st; readid<opt.read_ed; readid++) {
    int seq_len = strlen(readfile[readid]);
    if(seq_len>max_seq_len) max_seq_len = seq_len;
  }
  print_time("found max_seq_len");

  // candidate hypothesis
  vector<tr1::tuple<int, int, int> > hypothesis;
  generateHypothesis(opt.h_rate>0, hypothesis);
  
  print_time("generated hypothesis");
  // Confusion matrix.
  // Error loglikelihood matrix.
  double confMat[max_seq_len][4][4];
  double loglikelihood[max_seq_len][4][4];
  double*** confMat1 = new double**[max_seq_len];
  for (int i = 0; i < max_seq_len; ++i){
    confMat1[i] = new double*[4];
    for (int j = 0; j < 4; ++j){
      confMat1[i][j] = new double[4];
    }
  }
  initLoglikelihoodMat(opt, max_seq_len, confMat, loglikelihood, confMat1);
  print_time("init matrix");

  // Voting Mechanism.
  // Statitstics.
  vector<int> histogram(500, 0);
  set<unsigned int> hist_readset;

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



  // Zero out confusion matrix.
  zeroMat(max_seq_len, confMat);

  for(unsigned int readid = opt.read_st; readid<opt.read_ed; readid++) {
    //if (readid != 1)
    //	continue;
    print_time("begin collect votes for readid" + readid);
    string orig_seq = string(readfile[readid]);
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
      fout << corr_seq << endl;
      fqual << qual_seq << endl;
      continue;
    }

    // construct neighbors
    tr1::shared_ptr<map<unsigned int, NeighborInfo> > readNeighbors(new map<unsigned int, NeighborInfo>);
    for(size_t fiter=0; fiter<opt.inputFNames.size(); fiter++) {
      tr1::shared_ptr<map<unsigned int, NeighborInfo> > newNeighbors = neighborLoader[fiter]->get(readid);
      for(map<unsigned int, NeighborInfo>::iterator nn = newNeighbors->begin(); nn!=newNeighbors->end(); nn++) {
        map<unsigned int, NeighborInfo>::iterator conflict = readNeighbors->find(nn->first);
        int seq_len2 = strlen(readfile[nn->first]);
        if(conflict==readNeighbors->end() || // No conflict.
            ( nn->second!=conflict->second && nn->second.isBetter(conflict->second, seq_len, seq_len2) ))// Different alignment and conflict resolution.
          (*readNeighbors)[nn->first] = nn->second;
      }
    }
    print_time("constructed neighbors"); 
    vector<Brother> neighborsForClust;

    //find realNeigbor by Hammer
    for(map<unsigned int, NeighborInfo>::iterator neighbors = readNeighbors->begin();
        neighbors!=readNeighbors->end(); neighbors++) {
      const unsigned int& neighborId = neighbors->first;
      const char* neighbor_seq = readfile[neighborId];
      const int neighbor_seq_len = strlen(neighbor_seq);
      const int overlap = neighbors->second.get_overlap(seq_len, neighbor_seq_len);
      if(!neighbors->second.isNeighbor(seq_len, neighbor_seq_len, opt.K, opt.h, opt.e))
        continue;
      ostringstream new_brother;
      if (neighbors->second.get_offset() < 0){              
        for (int iter = -neighbors->second.get_offset(); iter < seq_len; ++iter){
          new_brother << neighbor_seq[iter];
        }  
        for (int iter = 0; iter < -neighbors->second.get_offset(); ++iter){
          new_brother << '_';
        }

      } else {
        for (int iter = 0; iter <neighbors->second.get_offset(); ++iter){
          new_brother << '_';
        }
        for (int iter = neighbors->second.get_offset(); iter < seq_len; ++iter){
          new_brother << neighbor_seq[iter - neighbors->second.get_offset()];
        }

      }
      neighborsForClust.push_back(Brother(neighborId, new_brother.str()));
    }
    print_time("make neighbors in nessesary format " + neighborsForClust.size());
    vector<Brother> realNeighbors = findRealBrothers(orig_seq, neighborsForClust, confMat1);
    print_time("found real neighbors");
    tr1::shared_ptr<map<unsigned int, NeighborInfo> > realNeighborsMap(new map<unsigned int, NeighborInfo>);
    for (size_t real_neighbor_id = 0; real_neighbor_id < realNeighbors.size(); ++real_neighbor_id){
      unsigned int id =(unsigned int) realNeighbors[real_neighbor_id].read_id;
      (*realNeighborsMap)[id] = (*readNeighbors)[id];
    }
    readNeighbors = realNeighborsMap;
    print_time("begin real voting");

    // Voting.
    // Collect votes.
    set<unsigned int> my_neighbors;
    for(map<unsigned int, NeighborInfo>::iterator neighbors = readNeighbors->begin();
        neighbors!=readNeighbors->end(); neighbors++) {

      const unsigned int& neighborId = neighbors->first;
      const char* neighbor_seq = readfile[neighborId];
      const int neighbor_seq_len = strlen(neighbor_seq);
      const int overlap = neighbors->second.get_overlap(seq_len, neighbor_seq_len);
      const double log_quality = 0;

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
      for(int offset=0; offset<overlap; offset++) {
        int calledb = baseToInt(overlap_seq[offset]);
        if(!isAnyBase(calledb)) {
          if(!orig) calledb = 3-calledb;
          Votes[st1+offset].push_back(VoteInfo(pos, calledb, log_quality, !orig));
          nVotes[st1+offset] += 1;
          // Original read gets to vote twice.
          if(readid==neighborId) {
            int nprior = max(1, (int)floor(sqrt(opt.cov)));
            if(opt.h_rate!=0) nprior = 0;
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
        for(vector<tr1::tuple<int,int,int> >::iterator H=hypothesis.begin(); H!=hypothesis.end(); H++) {
          int& b1 = tr1::get<0>(*H);
          int& b2 = tr1::get<1>(*H);
          int& b1b2 = tr1::get<2>(*H);
          int read_b1 = pVote->reverse_complement?3-b1:b1;
          int read_b2 = pVote->reverse_complement?3-b2:b2;

          if(b1>b2) {
            base_loglikelihood[b1b2] = -numeric_limits<float>::infinity();
            if(!pVote->prior)
              base_logquality[b1b2] = -numeric_limits<float>::infinity();			
          } else if(b1==b2)  {
            // Homozygous case.
            base_loglikelihood[b1b2] += loglikelihood[pVote->pos][pVote->base][read_b1] + pVote->log_quality;
            if(!pVote->prior)			
              base_logquality[b1b2] += loglikelihood[pVote->pos][pVote->base][read_b1] + pVote->log_quality;			
          } else  { // b1 < b2
            // Heterozygous case.
            base_loglikelihood[b1b2] += log(0.5*exp(loglikelihood[pVote->pos][pVote->base][read_b1]) + 0.5*exp(loglikelihood[pVote->pos][pVote->base][read_b2])) + pVote->log_quality;
            if(!pVote->prior)			
              base_logquality[b1b2] += log(0.5*exp(loglikelihood[pVote->pos][pVote->base][read_b1]) + 0.5*exp(loglikelihood[pVote->pos][pVote->base][read_b2])) + pVote->log_quality;
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
      if(total_prob>0 && !isAnyBase(orig_base))
        for(vector<tr1::tuple<int,int,int> >::iterator H=hypothesis.begin(); H!=hypothesis.end(); H++) {
          int& b1 = tr1::get<0>(*H);
          int& b2 = tr1::get<1>(*H);
          int& b1b2 = tr1::get<2>(*H);

          float prob = exp(base_logquality[b1b2]-base_logquality[tr1::get<2>(mse_base)]);
          if(readfile.isOrig(readid)) {
            confMat[i][b1][orig_base]+=0.5*prob/total_prob;
            confMat[i][b2][orig_base]+=0.5*prob/total_prob;    
          } else  {
            confMat[seq_len-i-1][3-b1][3-orig_base]+=0.5*prob/total_prob;
            confMat[seq_len-i-1][3-b2][3-orig_base]+=0.5*prob/total_prob;			    
          }
        }
        print_time("end real voting");
    }

    // Update histogram.
    if(opt.save_stats){
      size_t max_vote = *max_element(nVotes.begin(), nVotes.end());
      if(max_vote>histogram.size())
        histogram.resize(max_vote, 0);
      if(hist_readset.find(readid)==hist_readset.end()){
        // All neighbors will not vote in the future and only beginning of the read gets to vote to improve independence.
        hist_readset.insert(my_neighbors.begin(), my_neighbors.end());
        if(readfile.isOrig(readid))
          histogram[nVotes[0]]+=1;
        else
          histogram[nVotes[seq_len-1]]+=1;
      }
    }

    // Output corrected read.
    fout << corr_seq << endl;
    fqual << qual_seq << endl;
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
  ostringstream fconfMat;
  fconfMat <<"confMat.txt";
  ofstream fMat;
  fMat.open(fconfMat.str().c_str());
  for(int pos=0; pos<max_seq_len; pos++) {
    for(int b1=0; b1<4; b1++) {
      for(int b2=0; b2<4; b2++)
        fMat << confMat[pos][b1][b2] << ' ';
      fMat << endl;
    }
    fMat << endl;
  }
  fMat.close();

  return 0;
}
