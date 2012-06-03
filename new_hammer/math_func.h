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
#include "brother.h"
using namespace std;

double lMultinomial(const vector<Brother>& kmers);

double lBeta(vector<vector<Brother> >& groups);

int loaded_dice(vector<double> ps);

#endif /* MATH_FUNC_H_ */
