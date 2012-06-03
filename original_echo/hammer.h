#ifndef _HAMMER_H
#define _HAMMER_H

#include "brother.h"
#include <vector>
#include<string>

using namespace std;

vector<Brother> findRealBrothers(string read, vector<Brother>& brothers, double*** confMat);
#endif
