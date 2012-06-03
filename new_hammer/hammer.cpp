#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <utility>
#include <math.h>
#include <time.h>
#include <limits>
#include "DNASeq.hpp"
#include "math_func.h"
#include "collection_utils.h"
#include "distances.h"
#include "clusters.h"
#include "brother.h"

double EPS = 0.0000000000001;

using namespace std;

void print_array(vector<string> arr){
	cout<<"array size = "<<arr.size() << "\n";
	for (size_t i = 0; i < arr.size(); ++i){
		cout << arr[i]<<" ";
	}
	cout<<"\n\n";
}
double logLikelihood(Clusters& clusters, double*** confMat) {
	if (clusters.centers.size() == 0){
		return 0;
	}
	double logLikelihood = 0;
	if (clusters.centers.size() != 1){
		//logLikelihood += lgamma(clusters.centers.size());
		//logLikelihood += lBeta(clusters.groups);
	}
	for (size_t i = 0; i < clusters.centers.size(); ++i){
		for (size_t j = 0; j < clusters.groups[i].size(); ++j){
			logLikelihood += log(1- dist(clusters.centers[i], clusters.groups[i][j].seq, confMat))
					+ lMultinomial(clusters.groups[i]);
		}
	}
	return logLikelihood;
}

string prep_centroid(string centroid, vector<Brother>& reads, double*** confMat){
	if (centroid.find('_') != string::npos) {
		return centroid;
	}
	vector<double> probs_nonnorm;
	for (size_t i = 0; i < reads.size(); ++i){
		probs_nonnorm.push_back(1 - dist(centroid, reads[i].seq, confMat));
	}

	for (size_t i = 0; i < centroid.length(); ++i){
		if (centroid[i] != '_'){
			continue;
		}
		vector<double> probs_nonnorm_copy(probs_nonnorm);
		for (size_t read_idx = 0; read_idx < reads.size(); ++read_idx){
			if (reads[read_idx].seq[i] == '_'){
				probs_nonnorm_copy[read_idx] = 0;
			}
		}
		if (sum_vector(probs_nonnorm_copy) < EPS){
			continue;
		}
		int donor_idx = loaded_dice(normalize (probs_nonnorm_copy));
		centroid[i] = reads[donor_idx].seq[i];
	}
	return centroid;
}


vector<string> kpp(int k, vector<Brother>& reads, double*** confMat){
    srand(time(0));
    vector<string> centroids;
    string centroid = prep_centroid(reads[rand() % reads.size()].seq, reads, confMat);
    centroids.push_back(centroid);
    while ((int)centroids.size() < k ){
    	vector<double> sqrd_distances;
    	for (size_t l = 0; l < reads.size(); ++l){
    		double min_distance = min_dist(reads[l].seq, centroids, confMat).second;
    		sqrd_distances.push_back(min_distance * min_distance);
    	}
    	vector<double> probs = normalize(sqrd_distances);
    	int new_centroid_idx = loaded_dice(probs);
    	if (new_centroid_idx< 0 || new_centroid_idx > (int) reads.size()-1){
    		cerr<< "loaded_dice incorrect" + new_centroid_idx;
    	}
    	string new_centroid = prep_centroid(reads[new_centroid_idx].seq, reads, confMat);
    	bool was = false;
    	for (size_t j = 0; j < centroids.size(); ++j){
    		string oldcenter = centroids[j];
    		if (new_centroid == centroids[j]){
    			was = true;
    			break;
    		}
    	}
    	if (!was){
    		centroids.push_back(new_centroid);
    	}
    }
    return centroids;
}

string find_centroid(vector<Brother> reads, double*** confMat){
    string centroid;
    if (reads.size() == 0){
    	cout<<"error";
    	return "";
    }
    for (size_t i = 0; i < reads[0].seq.length(); ++i){
    	double probs[4] = {1.0, 1.0, 1.0, 1.0};
    	for (size_t j = 0; j < reads.size(); ++j){
    		string read = reads[j].seq;
    		char base = read[i];
    		if (base == '_'){
    			continue;
    		}
    		int base_idx = baseToInt(base);
    		if (base_idx < 4){
				for (size_t l = 0; l < 4; ++l){
					double observed = confMat[i][base_idx][l];
					probs[l] *= observed;
				}
    		}
    	}
    	int most_likely_idx = 0;
    	double max = 0;
    	for (int j = 0; j < 4; ++j){
    		if (probs[j] > max){
    			max = probs[j];
    			most_likely_idx = j;
    		}
    	}
    	centroid += intToBase(most_likely_idx);
    }
    return centroid;
}


Clusters* kmeans(int k, vector<Brother>& reads, double*** confMatrix){
	vector<string> centroids(kpp(k, reads, confMatrix));
	vector<string> prev_centroids = centroids;
	vector<vector<Brother> >  groups;
	bool first_pass = true;
	int iteration = 1;
	while ((first_pass || !vector_equals(centroids, prev_centroids)) && iteration < 1000){
		first_pass = false;
		prev_centroids = centroids;
		vector<int> assignments;
		for (size_t i = 0; i < reads.size(); ++i){
			assignments.push_back(min_dist(reads[i].seq, centroids, confMatrix).first);
		}
		groups.clear();
		for (size_t i = 0; i < centroids.size(); i++){
			vector<Brother> group;
			groups.push_back(group);
		}
		for (size_t i = 0; i < assignments.size(); i++){
			groups[assignments[i]].push_back(reads[i]);
		}
		centroids.clear();
		for (size_t i = 0; i < groups.size(); ++i){
			centroids.push_back(find_centroid(groups[i], confMatrix));
		}
		iteration++;
	}
	return new Clusters(centroids, groups);
}

vector<Brother> findRealBrothers(string read, vector<Brother>& brothers, double*** confMat){
	int num_clust = 1;
	double prevLog = -numeric_limits<double>::infinity( ) ;
	double log = -numeric_limits<double>::infinity( );
	Clusters* clusters;
	do{
		prevLog = log;
		clusters = kmeans(num_clust, brothers, confMat);
		log = logLikelihood(*clusters, confMat);
		num_clust++;
		cout << log << "\n";
	} while(prevLog < log);
	clusters = kmeans(num_clust - 2, brothers, confMat);
	int num = min_dist(read, clusters->centers, confMat).first;

	return clusters->groups[num];
}

int main(){
	//int k  = 2;
	vector<Brother> reads;
	string test_read = "GCGCTAC";
	Brother brother1(1, "___ATGC");
	reads.push_back(brother1);
	reads.push_back(Brother(1, "GCGCTGC"));
//	reads.push_back("GCGCTAC");
	reads.push_back(Brother(2, "GCGCT__"));
	reads.push_back(Brother(3, "TATATAT"));
	reads.push_back(Brother(4, "TCTAT__"));
	reads.push_back(Brother(5, "TATGTAT"));

	double** confMatBase = new double*[4];
	for (size_t i = 0; i < 4; ++i){
		confMatBase[i] = new double[4];
		for (int j = 0; j < 4; ++j)
			confMatBase[i][j] = 0.001;
		confMatBase[i][i] = 0.997;
	}
	double*** confMat = new double**[reads[0].seq.length()];
	for (size_t i = 0; i < reads[0].seq.length(); ++i){
		confMat[i] = confMatBase;
	}
//	vector<string> centroids = kmeans(k, reads, confMat)->centers;
//	for (size_t i = 0; i < centroids.size(); ++i){
//		cout << centroids[i]<<"\n";
//	}

	vector<Brother> brothers = findRealBrothers(test_read, reads, confMat);
	for (size_t i = 0; i < brothers.size(); ++i){
		cout << brothers[i].seq<<"\n";
	}
	return 0;
}
