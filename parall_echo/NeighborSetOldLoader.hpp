#ifndef __NEIGHBORSETOLD_H__
#define __NEIGHBORSETOLD_H__
#include "NeighborSet.hpp"
class NeighborSetOldLoader {
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
    NeighborSetOldLoader(char* NeighborFName) {
        cacheInit = true;
        cacheSt = 0;
        NeighborFile = fopen(NeighborFName, "rb");
    }

    ~NeighborSetOldLoader() {
        if(NeighborFile!=NULL)
            fclose(NeighborFile);
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
#endif
