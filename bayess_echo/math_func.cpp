/*
 * math_func.cpp
 *
 *  Created on: Jun 1, 2012
 *      Author: ira
 */

#include "math_func.h"
#include <cmath>
#include <map>
#include <time.h>
#include <stdlib.h>
#include "brother.h"

double lMultinomial(const vector<Brother>& kmers) {
  map<string, int> kmers_count;
  for (size_t i = 0; i < kmers.size(); ++i){
    if (kmers_count.find(kmers[i].seq) == kmers_count.end()){
      kmers_count.insert(pair<string, int>( kmers[i].seq, 1));
    } else{
      kmers_count[kmers[i].seq] += 1;
    }
  }
  double res = 0.0, sum = 0.0;
  for (map<string, int>::iterator it = kmers_count.begin(); it != kmers_count.end(); ++it){
    res += lgamma(it->second);
    sum += it->second;
  }
  return lgamma(sum) - res;
}

double lBeta(vector<vector<Brother> >& groups) {
  double res = 0.0, mult = 1.0;
  for (size_t i = 0; i < groups.size(); ++i) {
    res += lgamma(groups[i].size());
    mult *= groups[i].size();
  }
  return (res - lgamma(mult));
}

int loaded_dice(vector<double> ps){
  vector<double> cum_ps;
  cum_ps.push_back(ps[0]);
  for (size_t i = 1; i < ps.size(); ++i){
    cum_ps.push_back(cum_ps[i-1] + ps[i]);
  }
  double t = rand()/ double(RAND_MAX);
  for (size_t i = 0; i < ps.size(); ++i){
    if (t < cum_ps[i]){
      return i;
    }
  }
  return -1;
}
