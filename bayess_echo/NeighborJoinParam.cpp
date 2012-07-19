#include <iostream>
#include <stdexcept>
#include <vector>
#include <tr1/unordered_map>
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
#include "NeighborSet.hpp"

using namespace std;

class removeNeighborPredicate {
    public:
        removeNeighborPredicate(NeighborSet &Neighbors) : Neighbors(Neighbors) {};

        bool operator()(KmerOccurrence const &occ) {
            std::map<unsigned int, std::map<unsigned int, NeighborInfo> >::iterator iter_neighbor = Neighbors.getData().find(occ.readID);
            if(iter_neighbor==Neighbors.getData().end()) {
                return false;
            } else {
                return iter_neighbor->second.size()>NeighborSet::MaxNNeighbors;
            }
        }

        NeighborSet &Neighbors;
};

int main(int argc, char** argv) {
    // Initialize constants.
    Options opt(argc, argv);

    unsigned int read_st = opt.read_st;
    unsigned int read_ed = opt.read_ed;

    unsigned int ihash_st = opt.ihash_st;
    unsigned int ihash_ed = opt.ihash_ed;

    // Read data.
    // Check for files first.
    FILE *tmpfile = fopen(opt.inputFNames[0], "rb");
    if(tmpfile == NULL) {
        throw std::runtime_error("Could not open file 1.");
    } else {
        fclose(tmpfile);
    }

    tmpfile = fopen(opt.inputFNames[1], "rb");
    if(tmpfile == NULL) {
        throw std::runtime_error("Could not open file 2.");
    } else {
        fclose(tmpfile);
    }

    // Input files are hash files.
    HashMMAP kmer_mmap(opt.inputFNames[0], opt.inputFNames[1], ihash_st, ihash_ed);

    // Initialize new Neighbor Set (adjacency list).
    NeighborSet Neighbors(opt.read_st, opt.read_ed, opt.read_st, opt.read_ed);
    MMAPReads readfile(opt.readFName);

    // Construct Neighbor Sets
    for(HashMMAP::ConstIterator kmer_iter = kmer_mmap.begin(); kmer_iter!=kmer_mmap.end(); ++kmer_iter) {
        tr1::unordered_map<unsigned int, DNASeq> Reads;

        // Initialize neighbors.
        vector<KmerOccurrence> kmer_occurrences;
        vector<KmerOccurrence> kmer_occurrences2;

        for(HashMMAP::ConstReadIterator read=kmer_iter.read_begin(); read!=kmer_iter.read_end(); ++read) {
            // Add current read to hash from readIDs to read strings
            if(read.getReadID() >= read_st && read.getReadID() < read_ed) {
                Reads[read.getReadID()] = readfile[read.getReadID()];
                kmer_occurrences.push_back(KmerOccurrence(read.getReadID(), read.getPos()));
            }
            Reads[read.getReadID()] = readfile[read.getReadID()];
            kmer_occurrences2.push_back(KmerOccurrence(read.getReadID(), read.getPos()));
        }

        // Remove every read that has more than NeighborSet::MaxNNeighbors
        //kmer_occurrences.erase(remove_if(kmer_occurrences.begin(), kmer_occurrences.end(), [&Neighbors](const KmerOccurrence& occ) {
        //            auto iter_neighbor = Neighbors.getData().find(occ.readID);
        //            if(iter_neighbor==Neighbors.getData().end()) {
        //                return false;
        //            } else {
        //                return iter_neighbor->second.size()>NeighborSet::MaxNNeighbors;
        //            }
        //        }), kmer_occurrences.end());
        kmer_occurrences.erase(remove_if(kmer_occurrences.begin(), kmer_occurrences.end(), removeNeighborPredicate(Neighbors)), kmer_occurrences.end());

        // Compute neighbor set
        for(vector<KmerOccurrence>::iterator read1=kmer_occurrences.begin(); read1!=kmer_occurrences.end(); ++read1) {
            // Add neighbors from read1 to neighbor vertices
            const unsigned int readid1 = read1->readID;
            for(vector<KmerOccurrence>::iterator read2=kmer_occurrences2.begin(); read2!=kmer_occurrences2.end(); ++read2) {
                const unsigned int readid2 = read2->readID;
                if(readid2<readid1) continue;

                map<unsigned int, map<unsigned int, NeighborInfo> >::iterator iter_neighbor1 = Neighbors.getData().find(readid1);
                if(iter_neighbor1 != Neighbors.getData().end()) {
                    if(iter_neighbor1->second.size()>NeighborSet::MaxNNeighbors) break;
                }

                NeighborInfo old_ninfo, new_ninfo;
                new_ninfo.set_offset(read1->pos, read2->pos);

                map<unsigned int, NeighborInfo>::iterator iter_old_ninfo = Neighbors.getData()[readid1].find(readid2);
                bool seen = (iter_old_ninfo!=Neighbors.getData()[readid1].end());
                if( seen ) old_ninfo = iter_old_ninfo->second;
                if( seen && new_ninfo==old_ninfo) continue;

                const DNASeq& seq1 = Reads[readid1];
                const DNASeq& seq2 = Reads[readid2];

                bool is_neighbor = new_ninfo.isNeighbor(seq1, seq2, opt.K, opt.h, opt.e);
                bool is_better = seen ? is_neighbor&&new_ninfo.isBetter(old_ninfo, seq1.size(), seq2.size()) : true;

                if((!seen && is_neighbor) || (seen && is_better)) {
                    Neighbors.getData()[readid1][readid2] = NeighborInfo(read1->pos, read2->pos, new_ninfo.get_nerr());
                    if(readid2 >= read_st && readid2 < read_ed) {
                        Neighbors.getData()[readid2][readid1] = NeighborInfo(read2->pos, read1->pos, new_ninfo.get_nerr());
                    }
                }
                // In case the neighbor already exists, choose the one with minimum mismatch.
            }
        }
    }
    Neighbors.dump(opt.fpre, opt.fsuf);

    return 0;
}
