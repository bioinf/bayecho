#ifndef __MMAP_H__
#define __MMAP_H__

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdexcept>
#include <iostream>
#include <fstream>


class MMAP {
  int fdes;
  unsigned long long flen;
  void* mem;

  public:
  MMAP(char* fname) {
    struct stat buf;

    // Check for file existence.
    FILE *tmpfile;
    tmpfile = fopen(fname, "rb");
    if(tmpfile == NULL) {
      throw std::runtime_error("Could not open file for memory map.");
    } else {
      fclose(tmpfile);
    }

    // Construct memory map.
    fdes = open(fname, O_RDONLY);
    fstat(fdes, &buf);
    flen = buf.st_size;
    mem = mmap(0, flen, PROT_READ, MAP_SHARED | MAP_FILE, fdes, 0);
  }

  virtual ~MMAP() {
    munmap(mem, flen);
    close(fdes);
  }

  inline unsigned long long size() const {
    return flen;
  }

  inline const void* operator[](unsigned long long pos) const {
    if(pos<flen)
      return (void*)((char*)mem+pos);
    else {
      std::cerr << "Out of MMAP bounds." << std::endl;
      std::cerr << "pos: " << pos << std::endl;
      std::cerr << "flen: " << flen << std::endl;
      return NULL;
    }
  }    
};

#endif
