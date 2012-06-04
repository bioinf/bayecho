/*
 * distances.h
 *
 *  Created on: Jun 1, 2012
 *      Author: ira
 */

#ifndef DISTANCES_H_
#define DISTANCES_H_
#include <string>
#include <vector>
using namespace std;


double dist(string a, string b, double*** confMat);

pair<int, double> min_dist(string a, vector<string> centroids, double*** confMat);


#endif /* DISTANCES_H_ */
