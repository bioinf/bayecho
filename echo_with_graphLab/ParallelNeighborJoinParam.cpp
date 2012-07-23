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
#include <omp.h>
#include <tbb/concurrent_hash_map.h>

#include "util.hpp"
#include "DNASeq.hpp"
#include "MMAPReads.hpp"
#include "KmerHashMap.hpp"
#include "NeighborSet.hpp"

using namespace std;

typedef tbb::concurrent_hash_map<unsigned int, std::map<unsigned int, NeighborInfo> > neighbors_map;

class removeNeighborPredicate {
    public:
        removeNeighborPredicate(const neighbors_map &neighbors) : neighbors(neighbors) {};

        bool operator()(KmerOccurrence const &occ) {
            neighbors_map::const_accessor ac;
            if (neighbors.find(ac, occ.readID)) {
                return ac->second.size() > NeighborSet::MaxNNeighbors;
            }
            ac.release();
            return false;
        }

        const neighbors_map &neighbors;
};

int main(int argc, char** argv) {
    // Initialize constants.
    Options opt(argc, argv);
   
    std::cout << "Parallel Neighboring (param)" << std::endl;
    std::cout << "Parameters: " << std::endl;
    std::cout << "\t read_st: " << opt.read_st << std::endl;
    std::cout << "\t read_ed: " << opt.read_ed << std::endl;
    std::cout << "\t ihash_st: " << opt.ihash_st << std::endl;
    std::cout << "\t ihash_ed: " << opt.ihash_ed << std::endl;
    std::cout << "\t reads: " << opt.readFName << std::endl;
    std::cout << "\t K: " << opt.K << std::endl;
    std::cout << "\t h: " << opt.h << std::endl;
    std::cout << "\t e: " << opt.e << std::endl;
    std::cout << "\t fpre: " << opt.fpre << std::endl;
    std::cout << "\t fsuf: " << opt.fsuf << std::endl;
    std::cout << "\t input: " << opt.inputFNames[0] << std::endl;
    std::cout << "\t input: " << opt.inputFNames[1] << std::endl << std::endl;

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
    //NeighborSet Neighbors(opt.read_st, opt.read_ed, opt.read_st, opt.read_ed);
    neighbors_map neighbors;
    
    MMAPReads readfile(opt.readFName);
 
    // Construct Neighbor Sets
    unsigned int nKmer = kmer_mmap.getKmerNum();
    unsigned long cnt_kmer = 0;
    unsigned long cnt_nei = 0;
    #pragma omp parallel
    {
        #pragma omp single
        {
            std::cout << "Using " << omp_get_num_threads() << " threads" << std::endl;
        }

        #pragma omp for
        for (int kmer_index = 0; kmer_index < nKmer; kmer_index++) { //OMP requires kmer_index to be signed
            HashMMAP::ConstIterator kmer_iter = kmer_mmap.iterator(kmer_index);
            tr1::unordered_map<unsigned int, DNASeq> Reads;

            // Initialize neighbors.
            vector<KmerOccurrence> kmer_occurrences;
            vector<KmerOccurrence> kmer_occurrences2;

            for (HashMMAP::ConstReadIterator read=kmer_iter.read_begin(); read != kmer_iter.read_end(); ++read) {
                // Add current read to hash from readIDs to read strings
                if(read.getReadID() >= read_st && read.getReadID() < read_ed) {
                    Reads[read.getReadID()] = readfile[read.getReadID()];
                    kmer_occurrences.push_back(KmerOccurrence(read.getReadID(), read.getPos()));
                }
                Reads[read.getReadID()] = readfile[read.getReadID()];
                kmer_occurrences2.push_back(KmerOccurrence(read.getReadID(), read.getPos()));
            }

            kmer_occurrences.erase(remove_if(kmer_occurrences.begin(), kmer_occurrences.end(), removeNeighborPredicate(neighbors)), kmer_occurrences.end());

            // Compute neighbor set
            for (vector<KmerOccurrence>::iterator read1 = kmer_occurrences.begin(); read1 != kmer_occurrences.end(); ++read1) {
                // Add neighbors from read1 to neighbor vertices
                const unsigned int readid1 = read1->readID;
                for (vector<KmerOccurrence>::iterator read2 = kmer_occurrences2.begin(); read2 != kmer_occurrences2.end(); ++read2) {
                    const unsigned int readid2 = read2->readID;
                    if (readid2 < readid1) continue;

                    neighbors_map::const_accessor ac_neighbor1;
                    if (neighbors.find(ac_neighbor1, readid1)) {
                        if (ac_neighbor1->second.size() > NeighborSet::MaxNNeighbors) break;
                    } else {
                        neighbors.insert(ac_neighbor1, readid1);
                    }

                    NeighborInfo old_ninfo, new_ninfo;
                    new_ninfo.set_offset(read1->pos, read2->pos);

                    map<unsigned int, NeighborInfo>::const_iterator iter_old_ninfo = ac_neighbor1->second.find(readid2);
                    bool seen = (iter_old_ninfo != ac_neighbor1->second.end());
                    if (seen) {
                        old_ninfo = iter_old_ninfo->second;
                    }

                    if (seen && new_ninfo == old_ninfo) {
                        continue;
                    }

                    const DNASeq& seq1 = Reads[readid1];
                    const DNASeq& seq2 = Reads[readid2];

                    bool is_neighbor = new_ninfo.isNeighbor(seq1, seq2, opt.K, opt.h, opt.e);
                    bool is_better = seen ? is_neighbor && new_ninfo.isBetter(old_ninfo, seq1.size(), seq2.size()) : true;

                    ac_neighbor1.release();

                    if ((!seen && is_neighbor) || (seen && is_better)) {
                        neighbors_map::accessor w_ac_neighbor1;
                        neighbors.find(w_ac_neighbor1, readid1);
                        w_ac_neighbor1->second[readid2] = NeighborInfo(read1->pos, read2->pos, new_ninfo.get_nerr());
			w_ac_neighbor1.release();
                        if(readid2 >= read_st && readid2 < read_ed) {
                            neighbors_map::accessor w_ac_neighbor2;
                            if(!neighbors.find(w_ac_neighbor2, readid2)) {
                                neighbors.insert(w_ac_neighbor2, readid2);
                            }
                            w_ac_neighbor2->second[readid1] = NeighborInfo(read2->pos, read1->pos, new_ninfo.get_nerr());
                            w_ac_neighbor2.release();
                        }
                    }
                    // In case the neighbor already exists, choose the one with minimum mismatch.
                }
            }
        }
    } //pragma

    //dump
    NeighborSetDumper NSDumper(opt.fpre, opt.fsuf);

    std::vector<unsigned int> read_ids;
    for (neighbors_map::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
        read_ids.push_back(it->first);
    }
    std::sort(read_ids.begin(), read_ids.end());

    for(std::vector<unsigned int>::const_iterator it = read_ids.begin(); it != read_ids.end(); ++it) {
	unsigned int read_id = *it;
        neighbors_map::const_accessor ac;
	neighbors.find(ac, read_id);
        if(ac->second.size() > 0) {
            NSDumper.dump(read_id, ac->second);
        }
    }

    return 0;
}
