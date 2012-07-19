#include "edge.h"

Edge::Edge(int offset1, int nerr1){
  this->offset = offset1;
  this->nerr = nerr1;  
}

void Edge::add_to_message(unsigned int index, char ch){
  message.add_to_message(index, ch);
}

Message& Edge::get_message(){
  return message;
}

