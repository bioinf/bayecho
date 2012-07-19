/*
 * collection_utils.cpp
 *
 *  Created on: Jun 1, 2012
 *      Author: ira
 */

#include "collection_utils.h"
#include <iostream>
double sum_vector(vector<double> array){
  double sum = 0;
  for (size_t i = 0; i < array.size(); ++i){
    sum += array[i];
  }
  return sum;
}

vector<double>& normalize(vector<double>& array){
  double sum = sum_vector(array);
  for (size_t i = 0; i < array.size(); ++i){
    array[i] /= sum;
  }
  return array;
}


bool vector_equals(vector<string>& a, vector<string>& b){
  if (a.size() != b.size()){
    return false;
  }
  for (size_t i = 0; i < a.size(); ++i){
    if (a[i] != b[i]){
      return false;
    }
  }
  return true;
}
