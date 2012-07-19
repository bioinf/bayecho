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
#include <iostream>
#include "brother.h"

using namespace std;

struct Clusters{
  vector<string> centers;
  vector<vector<Brother> > groups;

  Clusters(vector<string> p_centers, vector< vector<Brother> > p_groups){
    centers = p_centers;
    groups = p_groups;
  }

  void print(){
    cerr << "count_clusters "<< centers.size() << "\n";
    for (int i = 0; i < centers.size(); ++i){
      cout<< "center\n" << centers[i] <<"\n";
      cout << "group_size = " << groups[i].size()<<"\n";
      for (int j = 0; j < groups[i].size(); ++j){
        cout << groups[i][j].seq<< "\n";
      }
      cout <<"\n\n";
    }
  }
};

#endif /* CLUSTERS_CPP_ */
