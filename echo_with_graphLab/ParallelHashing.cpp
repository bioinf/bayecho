#include <iostream>
#include <sstream>
#include <vector>
#include <cstdio>
#include <algorithm>

#include <tr1/unordered_map>
#include <tr1/memory>
#include <omp.h>
#include <tbb/concurrent_hash_map.h>

#include "util.hpp"
#include "DNASeq.hpp"
#include "MMAPReads.hpp"
#include "KmerHashMap.hpp"

struct kmer_hash_compare {
    static size_t hash(const Kmer& a) { return my_hash(a); }
    static bool equal(const Kmer& a, const Kmer& b) { return my_cmp(a, b); }
private:
    Kmer_hash my_hash;
    Kmer_cmp my_cmp;
};

typedef std::vector<KmerOccurrence> occ_list_type;
typedef tbb::concurrent_hash_map<Kmer, occ_list_type, kmer_hash_compare> kmers_map_type;
typedef tr1::unordered_map<Kmer, KmerOccurrence, Kmer_hash, Kmer_cmp> inner_map_type;

void dump_hash(const kmers_map_type&, const char*, const char*);

int main(int argc, char** argv) {
    Options opt(argc, argv);
    
    kmers_map_type kmers_map;
    
    MMAPReads readfile(opt.readFName);
    
    #pragma omp parallel
    {
        #pragma omp single
        {
            std::cout << "Using " << omp_get_num_threads() << " threads" << std::endl;
        }
        
        #pragma omp for
        for(int read_index = opt.read_st; read_index < opt.read_ed; read_index++) {
            bool orig = readfile.isOrig(read_index);
            inner_map_type occ_map;
            tr1::shared_ptr<DNASeq> cur_read(new DNASeq(readfile[read_index]));

            for (int pos = 0; pos <= cur_read->size() - opt.K; pos++) {
                Kmer kmer(cur_read, pos, opt.K);
                if (occ_map.find(kmer) == occ_map.end() || !orig) {
                    occ_map[kmer] = KmerOccurrence(read_index, pos);
                }
            }
            
            for(inner_map_type::iterator it = occ_map.begin(); it != occ_map.end(); it++) {
                kmers_map_type::const_accessor ac;
                if (kmers_map.find(ac, it->first)) {
                    if (ac->second.size() > KmerHashMap::MaxBinSize) {
                        continue;
                    }
                } else {
                    kmers_map.insert(ac, it->first);
                }
                ac.release();
                
                kmers_map_type::accessor acw;
                if (kmers_map.find(acw, it->first)) {
                    acw->second.push_back(it->second);
                } else {
                    std::cout << "ERROR: read is not found" << std::endl; 
                }
                acw.release();
            }
        }
    }
    
    dump_hash(kmers_map, opt.fpre, opt.fsuf);
    
    return 0;
}

struct index_less {
  inline bool operator()(const std::tr1::tuple<Kmer, occ_list_type* >& a,
                                const std::tr1::tuple<Kmer, occ_list_type* >& b) const {
      
    const Kmer& kmer1 = std::tr1::get<0>(a);
    const Kmer& kmer2 = std::tr1::get<0>(b);
      
    if (kmer1.hash < kmer2.hash) {
        return true;
    }
    
    if (kmer1.hash > kmer2.hash) {
        return false;
    }
    
    return (kmer1.kmer < kmer2.kmer);
  }    
};

void dump_hash(const kmers_map_type& kmers_map, const char* prefix, const char* suffix) {
    FILE* fout;
    std::ostringstream fname;
    fname << prefix << suffix << ".hash";
    fout = fopen(fname.str().c_str(), "wb");
    
    FILE* findexout;
    std::ostringstream findexname;
    findexname << prefix << suffix << ".index";
    findexout = fopen(findexname.str().c_str(), "wb");

    // compute index
    std::vector<std::tr1::tuple<Kmer, occ_list_type* > > index;
    index.push_back(std::tr1::make_tuple(Kmer("Z"), (occ_list_type*) 0));

    for (kmers_map_type::const_iterator it = kmers_map.begin(); it != kmers_map.end(); it++) {
        index.push_back(std::tr1::make_tuple(it->first, &(it->second)));
    }

    // Remove the (Kmer("Z"),0) tuple.
    index.erase(index.begin());
    
    std::sort(index.begin(), index.end(), index_less);

    // Write out index.
    // Format: kmer string, number of reads with kmer
    // Writes out:
    // num kmers in block, (kmer, num reads with kmer)+
    unsigned int nKmer = index.size();
    fwrite(&nKmer, sizeof(unsigned int), 1, findexout); // number of kmers
    for (unsigned int kid = 0; kid < nKmer; kid++) { // index of kmer
        const char* kmer = std::tr1::get<0>(index[kid]).getKmer().c_str();
        unsigned int nread = std::tr1::get<1>(index[kid])->size();
        fwrite(kmer, sizeof(char), strlen(kmer) + 1, findexout);
        fwrite(&nread, sizeof(unsigned int), 1, findexout);
    }

    // Write out reads.
    // Format: read ID, hit position
    // Writes out:
    // read id, pos within read
    unsigned int nKmer = index.size();
    for (unsigned int kid = 0; kid < nKmer; kid++) {
        const Kmer& kmer = std::tr1::get<0>(index[kid]);
        occ_list_type& occ_list = *std::tr1::get<1>(index[kid]);
        for (occ_list_type::const_iterator it = occ_list.begin(); it != occ_list.end(); it++) {
            const unsigned int& readid = it->readID;
            const unsigned int& pos = it->pos;
            fwrite(&readid, sizeof(unsigned int), 1, fout);
            fwrite(&pos, sizeof(unsigned int), 1, fout);		    
        }
    }

    fclose(fout[i]);
    fclose(findexout[i]);
}
