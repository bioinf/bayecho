/*
 * math_func.h
 *
 *  Created on: Jun 1, 2012
 *      Author: ira
 */

#ifndef _MY_MATH_FUNC_H_
#define _MY_MATH_FUNC_H_
#include <vector>
#include <string>
using namespace std;

double lMultinomial(const vector<string>& kmers);

double lBeta(vector<vector<string> >& groups);

int loaded_dice(vector<double> ps);

#endif /* MATH_FUNC_H_ */
