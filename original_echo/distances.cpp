/*
 * distances.cpp
 *
 *  Created on: Jun 1, 2012
 *      Author: ira
 */
#include "distances.h"
#include "DNASeq.hpp"

double EPS_CONF_MAT = 0.000001;

double dist(string a, string b, double*** confMat){
	double dist = 1;
	for (size_t i = 0; i < a.length(); ++i){
		int a_index = baseToInt(a[i]);
		int b_index = baseToInt(b[i]);
		if (a_index < 4 && b_index < 4){
			dist *= confMat[i][a_index][b_index] + EPS_CONF_MAT;
		}
	}
	return 1 - dist;
}

pair<int, double> min_dist(string a, vector<string> centroids, double*** confMat){
	double min = dist(a, centroids[0], confMat);
	int index_min = 0;
	for (size_t i = 1; i < centroids.size(); ++i){
		double cur_dist = dist(a, centroids[i], confMat);
		if (cur_dist < min){
			min = cur_dist;
			index_min = i;
		}
	}
	return make_pair(index_min, min);
}

