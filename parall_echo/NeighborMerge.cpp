#include <iostream>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <algorithm>
#include <tr1/memory>
#include <cmath>
#include <cstring>

#include "util.hpp"
#include "DNASeq.hpp"
#include "MMAPReads.hpp"
#include "NeighborSet.hpp"
#include "NeighborSetOldLoader.hpp"

using namespace std;

int main(int argc, char** argv) {
    // Initialize constants.
    Options opt(argc, argv);

    // Initialize Reads.
    MMAPReads readfile(opt.readFName);

    // Initialize output file.
    NeighborSetDumper NSDumper(opt.fpre, opt.fsuf);

    // Initialize NeighborSetLoaders
    vector<tr1::shared_ptr<NeighborSetOldLoader> > neighborLoader;
    for(size_t fiter=0; fiter<opt.inputFNames.size(); ++fiter)
        neighborLoader.push_back(tr1::shared_ptr<NeighborSetOldLoader>(new NeighborSetOldLoader(opt.inputFNames[fiter])));

    for(unsigned int readid = opt.read_st; readid<opt.read_ed; ++readid) {
        string orig_seq = string(readfile[readid]);
        int seq_len = orig_seq.size();

        // Construct neighbors
        tr1::shared_ptr<map<unsigned int, NeighborInfo> > readNeighbors(new map<unsigned int, NeighborInfo>);
        for(size_t fiter=0; fiter<opt.inputFNames.size(); fiter++) {
            tr1::shared_ptr<map<unsigned int, NeighborInfo> > newNeighbors = neighborLoader[fiter]->get(readid);
            for(map<unsigned int, NeighborInfo>::iterator nn = newNeighbors->begin(); nn!=newNeighbors->end(); nn++) {
                map<unsigned int, NeighborInfo>::iterator conflict = readNeighbors->find(nn->first);
                int seq_len2 = strlen(readfile[nn->first]);
                if(conflict==readNeighbors->end() || // No conflict
                    ( nn->second!=conflict->second && nn->second.isBetter(conflict->second, seq_len, seq_len2) )) { // Different alignment and conflict resolution.
                    (*readNeighbors)[nn->first] = nn->second;
                }
            }
        }
        NSDumper.dump(readid, *readNeighbors);
    }

    return 0;
}
