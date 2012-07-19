#include "message.h"

Message::Message(){
}

Message::Message(map<unsigned int, char> correct_read_chars1){
  this->correct_read_chars = correct_read_chars1;
}

Message::Message(const Message& message_to_copy){
  correct_read_chars = message_to_copy.correct_read_chars;
}

map<unsigned int, char>& Message::get_message(){
  return correct_read_chars;
}

void Message::add_to_message(unsigned int index, char ch){
  correct_read_chars.insert(pair< unsigned int, char> (index, ch));
}
