#ifndef _EDGE_H
#define _EDGE_H

#include <map>
using namespace std;

struct Edge{
  int offset;
  int nerr;  
  map<unsigned int, char> correct_read_chars; 
public:
  Edge(int offset1, int nerr1):offset(offset1), nerr(nerr1){};
  int get_offset(){return offset;};
  int get_nerr(){return nerr;};
  
};
#endif