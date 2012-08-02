#ifndef __MMAPREADS_H__
#define __MMAPREADS_H__

#include "MMAP.hpp"

class MMAPReads:public MMAP {
  unsigned long long nRead;
  unsigned long long* ReadIdx;
  char* ReadMem;
  public:
  MMAPReads(char* fname):MMAP(fname) {
    // Number of reads.
    nRead = *(unsigned long long*)MMAP::operator[](0);
    // Points to start of reads.
    ReadIdx = (unsigned long long*)MMAP::operator[](sizeof(unsigned long long));
    // Points to end of reads.
    ReadMem = (char*)MMAP::operator[](sizeof(unsigned long long)*(1+nRead));
  }

  virtual ~MMAPReads() {}

  inline unsigned long long size() const {
    return nRead;
  }

  inline const char* operator[](unsigned int readid) const {
    if(readid<nRead)
      return ReadMem + ReadIdx[readid] + 1;
    return NULL;
  }

  inline bool isOrig(unsigned int readid) const {
    if( *(ReadMem+ReadIdx[readid]) == 'Y')
      return true;
    else
      return false;
  }
};

#endif
