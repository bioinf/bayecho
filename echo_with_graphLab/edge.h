#ifndef _EDGE_H
#define _EDGE_H

#include <map>
#include "message.h"

using namespace std;

struct Edge{
  private:
    Message message;
    bool real_neighbor;
  public:
    int offset;
    int nerr;   
  public:
    Edge(int offset1, int nerr1);
    void not_neighbor();
    bool is_neighbor();
    int get_offset(){return offset;};
    int get_nerr(){return nerr;};
    void add_to_message(unsigned int index, char ch);
    Message& get_message();  
};
#endif
