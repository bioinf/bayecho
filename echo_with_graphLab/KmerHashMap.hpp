#ifndef __KMERHASHMAP_H__
#define __KMERHASHMAP_H__

#include <tr1/tuple>
#include <set>
#include <iostream>
#include <vector>
#include <sstream>
#include <memory>
#include <map>

#include <cstring>
#include <cstdio>
#include <cassert>

#include "DNASeq.hpp"
#include "MMAP.hpp"

struct KmerOccurrence {
    unsigned int readID, pos;

    KmerOccurrence(unsigned int readID, unsigned pos) : readID(readID),pos(pos) {}
    KmerOccurrence() {} 
};

class KmerHashMap {

    // Maps Kmers to (read_id, kmer_hit_pos) using Kmer_less comparator.
    //
    // *** Change hashing structure here ***
    // If you use a hashing structure that does not store the objects
    // in sorted order, then you will need to sort the kmers later
    // as noted below.
    std::map<Kmer, std::vector<KmerOccurrence>, Kmer_less> Data;

    std::vector<KmerOccurrence> dummy;
    unsigned int ihash_st, ihash_ed, nhash;

public:
    static const unsigned int MaxBinSize;
        
    KmerHashMap() : ihash_st(0),ihash_ed(1),nhash(1) {};
    KmerHashMap(unsigned int ihash_st, unsigned int ihash_ed, unsigned int nhash) : ihash_st(ihash_st), ihash_ed(ihash_ed), nhash(nhash) {};

    // Retrieve all reads that have the given kmer.
    std::vector<KmerOccurrence>& operator[](const Kmer& kmer) {
        if(dummy.size()>=MaxBinSize) dummy.clear();
        // Only accept if KmerHash falls within the right bucket
        // (which is determined below)
        // and the size of the bin is correct.
        // Otherwise, return a dummy vector.

        // 'batch' is which batch of blocks the kmer belongs to.
        // The hash table is divided into blocks.
        // The blocks are grouped into batches and processed together. 
        unsigned int batch = Kmer_hash()(kmer) % nhash;

        if(batch<ihash_st || batch>=ihash_ed) return dummy;
        std::vector<KmerOccurrence>& occ = Data[kmer];
        if(occ.size()>=MaxBinSize) return dummy;
        return occ;
    };

    void dump_binary(char* prefix, char* suffix) {
        // Number of blocks in a batch.
        assert(ihash_ed > ihash_st);
        const unsigned int nbatch = ihash_ed - ihash_st;
        FILE* fout[nbatch];
        FILE* findexout[nbatch];
        for(unsigned int i=0; i<nbatch; i++) {
            std::ostringstream fname;
            std::ostringstream findexname;
            // prefix is tmpDIR.
            // suffix is the starting read number.
            fname << prefix << suffix << ".hash";
            findexname << prefix << suffix << ".index";
            fout[i] = fopen(fname.str().c_str(), "wb");
            findexout[i] = fopen(findexname.str().c_str(), "wb");
        }

        // compute index
        std::vector<std::vector<std::tr1::tuple<Kmer, unsigned int> > > index(nbatch);
        for(unsigned int i=0; i<nbatch; i++)
            index[i].push_back(std::tr1::make_tuple(Kmer("Z"), 0));

        for(std::map<Kmer, std::vector<KmerOccurrence>, Kmer_less>::iterator iter=Data.begin(); iter!=Data.end(); iter++) {
            // ibatch is which block the kmer is in (misnomer, should be iblock instead)
            unsigned int ibatch = Kmer_hash()(iter->first) % nhash;
            assert(ibatch>=ihash_st && ibatch<ihash_ed);
            // fnum is the offset of the block within its batch.
            unsigned int fnum = ibatch - ihash_st;
            // (kmer, # reads that have kmer)
            index[fnum].push_back(std::tr1::make_tuple(iter->first, iter->second.size()));
        }

        // Remove the (Kmer("Z"),0) tuple.
        for(unsigned int i=0; i<nbatch; i++)
            index[i].erase(index[i].begin());

        // ***
        // NOTE: Might need to sort the kmers here. 
        // This is only necessary if we don't use a std::map
        // because iterating through std::map iterates in sorted order
        // ***

        // Write out index.
        // Format: kmer string, number of reads with kmer
        // nbatch is the number of blocks within a single batch.
        for(unsigned int i=0; i<nbatch; i++) {
            // Writes out:
            // num kmers in block, (kmer, num reads with kmer)+
            unsigned int nKmer = index[i].size();
            fwrite(&nKmer, sizeof(unsigned int), 1, findexout[i]); // number of kmers
            for(unsigned int kid=0; kid<nKmer; kid++) { // index of kmer
                const char* kmer = std::tr1::get<0>(index[i][kid]).getKmer().c_str();
                unsigned int nread = std::tr1::get<1>(index[i][kid]);
                fwrite(kmer, sizeof(char), strlen(kmer)+1, findexout[i]);
                fwrite(&nread, sizeof(unsigned int), 1, findexout[i]);
            }
        }

        // Write out reads.
        // Format: read ID, hit position
        for(unsigned int i=0; i<nbatch; i++) {
            // Writes out:
            // read id, pos within read
            unsigned int nKmer = index[i].size();
            for(unsigned int kid=0; kid<nKmer; kid++) {
                Kmer kmer = std::tr1::get<0>(index[i][kid]);
                std::vector<KmerOccurrence>& occur_list = Data[kmer];
                for(std::vector<KmerOccurrence>::iterator occ = occur_list.begin(); occ!=occur_list.end(); occ++) {
                    const unsigned int& readid = occ->readID;
                    const unsigned int& pos = occ->pos;
                    fwrite(&readid, sizeof(unsigned int), 1, fout[i]);
                    fwrite(&pos, sizeof(unsigned int), 1, fout[i]);		    
                }
            }
        }

        for(unsigned int i=0; i<nbatch; i++) {
            fclose(fout[i]);
            fclose(findexout[i]);
        }
    }
}; // end KmerHashMap

class HashMMAP {
    MMAP hash_mmap;
    MMAP index_mmap;

    unsigned int nKmer;

    unsigned long long fpos; // index pointer
    unsigned long long fpos2; // hash pointer

    unsigned long long fbegin; // start of hash block
    unsigned long long fend; // end of hash block

    unsigned long long f2begin; // start of hash block
    unsigned long long f2end; // end of hash block
    
    std::vector<unsigned long long> kmer_hash_offsets;

    public:
    class ConstReadIterator {
        MMAP const &hash_mmap;
        unsigned long long fpos2;

        public:
        ConstReadIterator(MMAP const& hash_mmap, unsigned long long fpos2) : hash_mmap(hash_mmap), fpos2(fpos2) {}

        inline unsigned int getReadID() const {
            if(hash_mmap[fpos2] == NULL) {
                std::cout << "ConstReadIterator:: getReadID() out of bounds." << std::endl;
            }
            return *(unsigned int*)hash_mmap[fpos2];
        }

        inline unsigned int getPos() const {
            if(hash_mmap[fpos2+sizeof(unsigned int)] == NULL) {
                std::cout << "ConstReadIterator: getPos() out of bounds." << std::endl;
            }
            return *(unsigned int*)hash_mmap[fpos2+sizeof(unsigned int)];
        }

        inline ConstReadIterator& operator++() {
            fpos2 += 2 * sizeof(unsigned int);
            return *this;
        }

        inline ConstReadIterator& operator++(int) {
            operator++();
            return *this;
        }

        inline bool operator==(ConstReadIterator const &a) const {
            return fpos2 == a.fpos2;
        }

        inline bool operator!=(ConstReadIterator const &a) const {
            return !(fpos2 == a.fpos2);
        }
    };

    class ConstIterator {
        MMAP const &hash_mmap;
        MMAP const &index_mmap; 

        unsigned long long fpos;
        unsigned long long fpos2;

        public:
        ConstIterator(MMAP const& hash_mmap, MMAP const& index_mmap, unsigned long long fpos, unsigned long long fpos2) : hash_mmap(hash_mmap), index_mmap(index_mmap), fpos(fpos), fpos2(fpos2) {}

        inline Kmer getKmer() const {
            if(index_mmap[fpos] == NULL) {
                std::cout << "getKmer() out of bounds." << std::endl;
            }
            return Kmer((char*)index_mmap[fpos]);
        }

        inline unsigned int getOccur() const {
            if(index_mmap[fpos] == NULL) {
                std::cout << "getOccur() out of bounds." << std::endl;
            }
            return *(unsigned int*)index_mmap[strlen((char*)index_mmap[fpos]) + 1 + fpos];
        }

        inline ConstReadIterator read_begin() const {
            return ConstReadIterator(hash_mmap, fpos2);
        }

        inline ConstReadIterator read_end() const {
            return ConstReadIterator(hash_mmap, fpos2 + 2*getOccur() * sizeof(unsigned int));
        }

        inline ConstIterator& operator++() {
            int offset = strlen((char*)index_mmap[fpos]) + 1;
            if(index_mmap[offset + fpos] == NULL) {
                std::cout << "operator++() out of bounds." << std::endl;
            }
            fpos2 +=  2 * (*(unsigned int*)index_mmap[offset + fpos]) * sizeof(unsigned int);
            fpos += offset + sizeof(unsigned int);
            return *this;
        }

        inline ConstIterator& operator++(int) {
            operator++();
            return *this;
        }

        // Equality depends only on where the index pointer (fpos) points to.
        inline bool operator==(ConstIterator const& a) const {
            return fpos == a.fpos;
        }

        inline bool operator!=(ConstIterator const& a) const {
            return !(fpos == a.fpos);
        }
    };

    HashMMAP(char* fname, char* findex) : hash_mmap(fname), index_mmap(findex), fpos(0), fpos2(0) {
        nKmer = *(unsigned int*)index_mmap[0];
        fpos = sizeof(unsigned int);
        unsigned int len_kmer = strlen((char*)index_mmap[sizeof(unsigned int)])+1;

        fbegin = sizeof(unsigned int);
        fend = index_mmap.size();
        f2begin = 0;
        f2end = hash_mmap.size();
        
        kmer_hash_offsets.resize(nKmer, 0);
        unsigned long long tmpfpos = sizeof(unsigned int);
        unsigned int cur_kmer_index = 0;
        unsigned long long curpos = 0;
        while (tmpfpos < fend) {
            kmer_hash_offsets[cur_kmer_index] = curpos;
            tmpfpos += len_kmer;
            unsigned int noccur = *(unsigned int*)index_mmap[tmpfpos];
            curpos += noccur * 2 * sizeof(unsigned int);
            tmpfpos += sizeof(unsigned int);
            cur_kmer_index++;
        }
        if (cur_kmer_index != nKmer) {
            std::cout << "ERROR: kmer_hash_offsets was not initialized correctly!" << std::endl;
        }
        assert(cur_kmer_index == nKmer);        
    }

    HashMMAP(char* fname, char* findex, unsigned int ihash_st, unsigned int ihash_ed) : hash_mmap(fname), index_mmap(findex), fpos(0), fpos2(0), fbegin(0), fend(0), f2begin(0), f2end(0) {
        nKmer = *(unsigned int*)index_mmap[0];
        
        kmer_hash_offsets.resize(nKmer, 0);
        
        // I assume all kmers are the same length so just get the
        // length of first kmer.
        unsigned int len_kmer = strlen((char*)index_mmap[sizeof(unsigned int)])+1;
        fpos = fbegin = (len_kmer+sizeof(unsigned int)) * ihash_st + sizeof(unsigned int);
        fend = (len_kmer+sizeof(unsigned int)) * ihash_ed + sizeof(unsigned int);
        assert(fend <= index_mmap.size());

        // Set fpos2 to the correct position.
        // Start from beginning of the index.
        unsigned long long tmpfpos = sizeof(unsigned int);
        unsigned int cur_kmer_index = 0;
        fpos2 = 0;
        while(tmpfpos < fbegin) {
            kmer_hash_offsets[cur_kmer_index] = fpos2;
            tmpfpos += len_kmer; // Skip the kmer string.
            // Read the number of occurrences.
            unsigned int noccur = *(unsigned int*)index_mmap[tmpfpos]; 
            fpos2 += noccur * 2 * sizeof(unsigned int);
            tmpfpos += sizeof(unsigned int);
            cur_kmer_index++;
        }
        assert(tmpfpos == fbegin);
        assert(cur_kmer_index == ihash_st);

        // Set f2begin
        f2begin = fpos2;

        // Find end position.
        // *Continue* from tmpfpos, which was set above.
        // (Do not change ordering of these while loops.)
        f2end = fpos2;
        while(tmpfpos < fend) {
            kmer_hash_offsets[cur_kmer_index] = f2end;
            tmpfpos += len_kmer; // Skip the kmer string.
            // Read the number of occurrences.
            unsigned int noccur = *(unsigned int*)index_mmap[tmpfpos];
            f2end += noccur * 2 * sizeof(unsigned int);
            tmpfpos += sizeof(unsigned int);
            cur_kmer_index++;
        }
        assert(f2end <= hash_mmap.size());
        if (cur_kmer_index != ihash_ed) {
            std::cout << "ERROR: kmer_hash_offsets was not initialized correctly!" << std::endl;
        }
        assert(cur_kmer_index == ihash_ed);
    }

    inline ConstIterator begin() const {
        return ConstIterator(hash_mmap, index_mmap, fbegin, f2begin);
    }
    
    inline ConstIterator iterator(unsigned int kmer_index) const {
        unsigned int len_kmer = strlen((char*)index_mmap[sizeof(unsigned int)])+1;
        unsigned long long findexpos = (len_kmer + sizeof(unsigned int)) * kmer_index + sizeof(unsigned int);
        unsigned long long fhashpos = kmer_hash_offsets[kmer_index];
        
        return ConstIterator(hash_mmap, index_mmap, findexpos, fhashpos);
    }

    inline ConstIterator end() const {
        return ConstIterator(hash_mmap, index_mmap, fend, f2end);
    }
    
    inline unsigned int getKmerNum() const {
        return nKmer;
    }
};

#endif
