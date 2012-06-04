#ifndef __NEIGHBORSET_H__
#define __NEIGHBORSET_H__

#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <tr1/memory>

#include <cstdio>
#include <cstring>
#include <cassert>

#include "MMAPReads.hpp"
#include "DNASeq.hpp"

class NeighborInfo {
    public:
    typedef char index_t;

    private:
    index_t offset, nerr;

    public:
    NeighborInfo() {
        offset = 0;
        nerr = -1;
    }

    NeighborInfo(int offset, int nerr) : offset(offset), nerr(nerr) {}

    NeighborInfo(int pos1, int pos2, int nerr) : offset(pos1-pos2), nerr(nerr) {}

    inline void set_offset(int pos1, int pos2) {
        offset = pos1-pos2;
    }

    inline int get_st1() const {
        return static_cast<int>(offset < 0 ? 0 : offset);
    }

    inline int get_st2() const {
        return static_cast<int>(offset >= 0 ? 0 : -offset);
    }

    inline int get_offset() const {
        return static_cast<int>(offset);
    }

    inline int get_nerr() const {
        return static_cast<int>(nerr);
    }

    inline int get_overlap(int seqlen1, int seqlen2) const {
        int len1 = seqlen1 - get_st1();
        assert(len1 > 0);
        int len2 = seqlen2 - get_st2();
        assert(len2 > 0);

        int overlap = len1 <= len2 ? len1 : len2;
        return overlap;
    }

    bool isNeighbor(const DNASeq& read1, const DNASeq& read2, int k, int h, float e);

    inline bool isNeighbor(int seqlen1, int seqlen2, int k, int h, float e) const {
        assert(nerr!=-1);
        int overlap = get_overlap(seqlen1, seqlen2);
        assert(overlap >= 0);
        if(overlap<h) return false;
        if(1.0*nerr/overlap>e) return false;
        return true;
    }

    inline bool isBetter(const NeighborInfo& info2, int seqlen1, int seqlen2) const {
        assert(nerr!=-1);
        const NeighborInfo& info1 = *this;

        return (1.0*info1.get_nerr()/info1.get_overlap(seqlen1, seqlen2)<
                1.0*info2.get_nerr()/info2.get_overlap(seqlen1, seqlen2));
    }

    inline bool operator==(const NeighborInfo& b) {
        return offset==b.offset;
    };

    inline bool operator!=(const NeighborInfo& b) {
        return offset!=b.offset;
    };
};

class NeighborSetDumper {
    FILE* fout;

    NeighborSetDumper(); // No default constructor.

    public:
    NeighborSetDumper(const char* prefix, const char* suffix) {
        std::ostringstream fname;
        fname << prefix << "neighbors_" << suffix << ".list";
        fout = fopen(fname.str().c_str(), "wb");
    };

    ~NeighborSetDumper() {
        if(fout!=NULL) fclose(fout);
    };

    inline void dump(unsigned int read_id, const std::map<unsigned int, NeighborInfo>& neighbors) {
        if(neighbors.size()==0) return;
        // Write out current read id.
        fwrite(&read_id, sizeof(unsigned int), 1, fout);

        // Write out number of neighbors in neighbor list.
        unsigned int neighbors_size = neighbors.size();
        fwrite(&neighbors_size, sizeof(unsigned int), 1, fout);

        // Write out all neighbors of current read id
        for(std::map<unsigned int, NeighborInfo>::const_iterator info=neighbors.begin(); info!=neighbors.end(); info++) {
            const unsigned int &neighborID = info->first;
            const NeighborInfo& ninfo = info->second;

            NeighborInfo::index_t offset = ninfo.get_offset();
            NeighborInfo::index_t nerr = ninfo.get_nerr();

            // Write neighbor id.
            fwrite(&neighborID, sizeof(unsigned int), 1, fout);
            // Write offset.
            fwrite(&offset, sizeof(NeighborInfo::index_t), 1, fout);
            // Write number of errors.
            fwrite(&nerr, sizeof(NeighborInfo::index_t), 1, fout); 
        }
    }
};

class NeighborSetLoader {
    FILE* NeighborFile;
    std::map<unsigned int, std::tr1::shared_ptr<std::map<unsigned int, NeighborInfo> > > Neighbors;
    unsigned int cacheSt;
    bool cacheInit;

    static const unsigned int MaxSeqLen;
    static const unsigned int CacheSize;

    void load(unsigned int st, unsigned int ed) {
        if(NeighborFile==NULL) return;

        // Erase obsolete data.
        if( Neighbors.begin()->first < st )
            Neighbors.erase(Neighbors.begin(), Neighbors.lower_bound(st));

        // Load new data.
        unsigned int curID;
        unsigned int numNeighbors;
        unsigned int neighborId;
        NeighborInfo::index_t offset, nerr;

        // fread returns 0 when it can't read any more bytes.
        while(fread(&curID, sizeof(unsigned int), 1, NeighborFile)==1) {
            // Read number of neighbors.
            fread(&numNeighbors, sizeof(unsigned int), 1, NeighborFile);

            std::tr1::shared_ptr<std::map<unsigned int, NeighborInfo> > curNeighbors(new std::map<unsigned int, NeighborInfo>);

            for(unsigned int ii=0; ii<numNeighbors; ++ii) {
                // Read in neighbor id
                fread(&neighborId, sizeof(unsigned int), 1, NeighborFile);
                // Read in offset.
                fread(&offset, sizeof(NeighborInfo::index_t), 1, NeighborFile);
                // Read in number of errors.
                fread(&nerr, sizeof(NeighborInfo::index_t), 1, NeighborFile);

                curNeighbors->insert(std::make_pair(neighborId, NeighborInfo(offset, nerr)));
            }

            if(curID>=st)
                Neighbors[curID] = curNeighbors;
            if(curID>=ed)
                break;
        }
    }    

    public:
    NeighborSetLoader(char* NeighborFName) {
        cacheInit = true;
        cacheSt = 0;
        NeighborFile = fopen(NeighborFName, "rb");
    }

    ~NeighborSetLoader() {
        if(NeighborFile!=NULL)
            fclose(NeighborFile);
    }

    std::map<unsigned int, std::tr1::shared_ptr<std::map<unsigned int, NeighborInfo> > > get_map(int min_read_id, int max_read_id){
      load(min_read_id, max_read_id);
      return Neighbors;
    }

    inline std::tr1::shared_ptr<std::map<unsigned int, NeighborInfo> > get(unsigned int readID) {
        if(NeighborFile==NULL) return std::tr1::shared_ptr<std::map<unsigned int, NeighborInfo> >(new std::map<unsigned int, NeighborInfo>);

        if( cacheInit || readID >= cacheSt + CacheSize ) {
            while(readID > cacheSt + CacheSize) {
                if(cacheInit) {
                    cacheSt = 0;
                    cacheInit = false;
                } else {
                    cacheSt = cacheSt + CacheSize;
                }
            }
            load(cacheSt, cacheSt+CacheSize);
        }

        std::map<unsigned int, std::tr1::shared_ptr<std::map<unsigned int, NeighborInfo> > >::const_iterator it = Neighbors.find(readID);
        if(it != Neighbors.end())
            return it->second;
        else
            return std::tr1::shared_ptr<std::map<unsigned int, NeighborInfo> >(new std::map<unsigned int, NeighborInfo>);
    }
};

// NeighborSet is an adjacency list, where the std::map of source vertices stores a std::map of adjacent (destination) vertices. Attached to each adjacent vertex is the number of errors in the overlap and the position where the two reads overlap. 
class NeighborSet {
    unsigned int st1, ed1, range1;
    unsigned int st2, ed2, range2;
    std::map<unsigned int, std::map<unsigned int, NeighborInfo> > Data;
    std::map<unsigned int, NeighborInfo> dummy;

    NeighborSet() {};

    public:
    static const unsigned int MaxNNeighbors;

    NeighborSet(unsigned int st1, unsigned int ed1, unsigned int st2, unsigned int ed2) : st1(st1), ed1(ed1), st2(st2), ed2(ed2) {
        assert(st1 <= st2);
        assert(ed1 <= ed2);
        range1 = ed1-st1;
        range2 = ed2-st2;
        if(ed1-st2>0)
            range2 -= ed1-st2;
        assert(range1+range2 >= 0);
    }    

    void dump(char* prefix, char* suffix) const  {
        NeighborSetDumper NSDumper(prefix, suffix);

        for(std::map<unsigned int, std::map<unsigned int, NeighborInfo> >::const_iterator iter = Data.begin(); iter != Data.end(); ++iter) {
            if(iter->second.size() > 0) {
                unsigned int read_id = iter->first;
                NSDumper.dump(read_id, iter->second);
            }
        }
    }

    std::map<unsigned int, std::map<unsigned int, NeighborInfo> >& getData() {
        return Data;
    }
};

#endif
