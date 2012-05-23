#include <iostream>
#include <vector>
#include <tr1/unordered_map>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <utility>
#include <algorithm>
#include <tr1/memory>

#include "util.hpp"
#include "DNASeq.hpp"
#include "MMAPReads.hpp"
#include "KmerHashMap.hpp"

using namespace std;

int main(int argc, char** argv) {
    // Initialize constants.
    Options opt(argc, argv);

    // Input reads with reverse complements.
    KmerHashMap KmerBin(opt.ihash_st, opt.ihash_ed, opt.nhash);

    // Read data.
    MMAPReads readfile(opt.readFName);

    // Construct hash bin.
    for(unsigned int i=opt.read_st; i<opt.read_ed; i++) {
        bool orig = readfile.isOrig(i);
        tr1::unordered_map<Kmer, KmerOccurrence, Kmer_hash, Kmer_cmp> myKmer;
        tr1::shared_ptr<DNASeq> curRead(new DNASeq(readfile[i]));

        for(int pos=0; pos<=curRead->size()-opt.K; pos++) {
            Kmer kmer(curRead, pos, opt.K);
            if(myKmer.find(kmer)==myKmer.end()) {
                // new kmer -> insert
                myKmer[kmer]=KmerOccurrence(i, pos);
            } else  {
                // existing kmer -> keep the one that is closest to the beginning of the read
                if(!orig) myKmer[kmer] = KmerOccurrence(i, pos);
            }
        }

        // Insert kmer back to KmerBin
        for(tr1::unordered_map<Kmer, KmerOccurrence, Kmer_hash, Kmer_cmp>::iterator mykmer = myKmer.begin();
                mykmer != myKmer.end(); mykmer++) {
            KmerBin[mykmer->first].push_back(mykmer->second);
        }
    }

    // fpre is tmpDIR.
    // fsuf is the read start ID for the hash file.
    KmerBin.dump_binary(opt.fpre, opt.fsuf);

    return 0;
}
