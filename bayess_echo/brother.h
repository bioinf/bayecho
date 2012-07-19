/*
 * brother.h
 *
 *  Created on: Jun 3, 2012
 *      Author: ira
 */

#ifndef BROTHER_H_
#define BROTHER_H_
#include <string>

using namespace std;

struct Brother{
  int read_id;
  string seq;
  Brother(int m_read_id, string m_seq){
    read_id = m_read_id;
    seq = m_seq;
  }

  Brother(const Brother& brother){
    read_id = brother.read_id;
    seq = brother.seq;
  }
};

#endif /* BROTHER_H_ */
