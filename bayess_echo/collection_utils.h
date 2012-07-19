/*
 * collection_utils.h
 *
 *  Created on: Jun 1, 2012
 *      Author: ira
 */

#ifndef COLLECTION_UTILS_H_
#define COLLECTION_UTILS_H_
#include <vector>
#include <string>
using namespace std;

double sum_vector(vector<double> array);

vector<double>& normalize(vector<double>& array);

bool vector_equals(vector<string>& a, vector<string>& b);

#endif /* COLLECTION_UTILS_H_ */
