#include <cstring>

#include "DNASeq.hpp"

using namespace std;

tr1::unordered_map<char, char> DNASeq::ComplementBase = DNASeq::ComplementBaseInitializer();
