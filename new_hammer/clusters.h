/*
 * clusters.cpp
 *
 *  Created on: Jun 1, 2012
 *      Author: ira
 */

#ifndef CLUSTERS_CPP_
#define CLUSTERS_CPP_
#include <vector>
#include <string>
#include "brother.h"

using namespace std;

struct Clusters{
	vector<string> centers;
	vector<vector<Brother> > groups;

	Clusters(vector<string> p_centers, vector< vector<Brother> > p_groups){
		centers = p_centers;
		groups = p_groups;
	}
};

#endif /* CLUSTERS_CPP_ */
