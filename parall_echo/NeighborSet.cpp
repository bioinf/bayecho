#include "NeighborSet.hpp"
#include "DNASeq.hpp"

using namespace std;

const unsigned int NeighborSet::MaxNNeighbors = 100000;
const unsigned int NeighborSetLoader::MaxSeqLen = 500000;
const unsigned int NeighborSetLoader::CacheSize = 500000;

bool NeighborInfo::isNeighbor(const DNASeq& read1, const DNASeq& read2, int k, int h, float e) {
  int st1 = get_st1();
  int st2 = get_st2();
  int overlap = get_overlap(read1.size(), read2.size());

  if(overlap<h)
    return false;

  nerr = 0;
  const char* seq1 = read1.getSeq();
  const char* seq2 = read2.getSeq();

  for(int i=0; i<overlap; i++)
    if(seq1[st1+i]!=seq2[st2+i])
      nerr++;
  if(nerr>overlap*e)
    return false;
  return true;
}
