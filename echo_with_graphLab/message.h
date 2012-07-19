#ifndef _Message_H_
#define _Message_H_

#include <map>

using namespace std;

class Message{
  private:
    map<unsigned int, char> correct_read_chars;
  public:
    Message();
    Message(map<unsigned int, char> correct_read_chars);
    Message(const Message& message);
    map<unsigned int, char>& get_message();
    void add_to_message(unsigned int index, char ch); 
};

#endif
