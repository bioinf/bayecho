#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <utility>
#include <algorithm>
#include <memory>

#include "util.hpp"
#include "DNASeq.hpp"
#include "MMAPReads.hpp"
#include "KmerHashMap.hpp"

using namespace std;

int main(int argc, char** argv) {
  // Initialize input arguments.
  Options opt(argc, argv);

  assert(opt.inputFNames.size()%2 == 0);

  // Open hash files and index files
  std::vector<tr1::shared_ptr<MMAP> > hash_files;
  std::vector<tr1::shared_ptr<MMAP> > index_files;
  for(vector<char*>::iterator iter = opt.inputFNames.begin(); iter != opt.inputFNames.end(); ++iter) {
    hash_files.push_back(tr1::shared_ptr<MMAP>(new MMAP(*iter)));
    ++iter;
    index_files.push_back(tr1::shared_ptr<MMAP>(new MMAP(*iter)));
  }

  // Output hash table file.
  FILE *fout = fopen((string(opt.fpre) + opt.fsuf + ".hash").c_str(), "wb");

  // Output index file.
  FILE *findexout = fopen((string(opt.fpre) + opt.fsuf + ".index").c_str(), "wb");

  unsigned int nFiles = opt.inputFNames.size()/2;
  vector<unsigned long long> fpos(nFiles, 0); // points to index
  vector<unsigned long long> fpos2(nFiles, 0); // points to hash table

  vector<unsigned int> num_kmers_read(nFiles, 0);
  vector<unsigned int> nKmer(nFiles, 0);

  unsigned int num_tot_kmers = 0;

  for(unsigned int hash_id=0; hash_id<nFiles; ++hash_id) {
    fpos[hash_id] = 0;
    nKmer[hash_id] = *(unsigned int*)(*index_files[hash_id])[0];
    fpos[hash_id] += sizeof(unsigned int); // Advance file pointer.
    fpos2[hash_id] = 0;
  }

  // Leave a space for the total number of kmers at the beginning of index.
  unsigned int tmpint = 0;
  fwrite(&tmpint, sizeof(unsigned int), 1, findexout);

  // Iterate through all hash tables.
  bool all_end_reach = false;
  while(!all_end_reach) {
    // Check to see which hash files we've reached the end of.
    vector<unsigned int> good_hashes;
    for(unsigned int hash_id=0; hash_id<nFiles; ++hash_id) {
      if(num_kmers_read[hash_id] < nKmer[hash_id]) {
        good_hashes.push_back(hash_id);
      }
    }
    if(good_hashes.size() == 0) {
      all_end_reach = true;
      continue;
    }

    // Find "minimum" kmer among the hash files.
    // Initialize min_kmer to kmer in first good hash file.
    Kmer min_kmer((char*)(*index_files[*good_hashes.begin()])[fpos[*good_hashes.begin()]]);

    for(vector<unsigned int>::iterator iter=good_hashes.begin(); iter!=good_hashes.end(); ++iter) {
      Kmer kmer((char*)(*index_files[*iter])[fpos[*iter]]);
      if(Kmer_less()(kmer, min_kmer))
        min_kmer = string((char*)(*index_files[*iter])[fpos[*iter]]);
    }

    // Update total number of kmers.
    ++num_tot_kmers;

    // Find which hash tables have the minimum kmer
    // and adjust their fpos, etc. appropriately.
    vector<unsigned int> min_hashes;
    vector<unsigned int> nreads(nFiles);
    unsigned int tot_noccur = 0;
    for(vector<unsigned int>::iterator iter=good_hashes.begin(); iter != good_hashes.end(); ++iter) {
      Kmer kmer((char*)(*index_files[*iter])[fpos[*iter]]);
      if(Kmer_cmp()(kmer,min_kmer)) {
        min_hashes.push_back(*iter);
        Kmer kmer((char*)(*index_files[*iter])[fpos[*iter]]);
        fpos[*iter] += min_kmer.getKmer().size()+1; // Advance file ptr.

        // Get number of reads for the kmer.
        nreads[*iter]=*(unsigned int*)(*index_files[*iter])[fpos[*iter]];
        tot_noccur += nreads.at(*iter); 
        fpos[*iter] += sizeof(unsigned int); // Advance file ptr.
        ++num_kmers_read[*iter];
      }
    }

    const char* kmer_string = min_kmer.getKmer().c_str();
    // Write out index.
    //    Write kmer string.
    fwrite(kmer_string, sizeof(char), strlen(kmer_string)+1, findexout);
    //    Write number of reads with kmer.
    fwrite(&tot_noccur, sizeof(unsigned int), 1, findexout);

    // write out reads and positions
    for(vector<unsigned int>::iterator iter=min_hashes.begin(); iter != min_hashes.end(); ++iter) {
      for(unsigned int occ_count=0; occ_count<nreads[*iter]; ++occ_count) {
        unsigned int read_id = *(unsigned int*)(*hash_files[*iter])[fpos2[*iter]];
        fpos2[*iter] += sizeof(unsigned int);
        unsigned int pos = *(unsigned int*)(*hash_files[*iter])[fpos2[*iter]];
        fpos2[*iter] += sizeof(unsigned int);
        fwrite(&read_id, sizeof(unsigned int), 1, fout);
        fwrite(&pos, sizeof(unsigned int), 1, fout);
      }
    }
  }

  // Write out total number of kmers to beginning of file.
  fseek(findexout, 0L, SEEK_SET);
  fwrite(&num_tot_kmers, sizeof(unsigned int), 1, findexout);

  // Close files.
  fclose(findexout);
  fclose(fout);
}
