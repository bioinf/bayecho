#ifndef UPDATE_FUNCTION_
#define UPDATE_FUNCTION_
#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>
#include "node.h"
#include "edge.h"
typedef graphlab::graph<MyNode, Edge> graph_type;
typedef graphlab::types<graph_type> gl_types;

struct VoteInfo {
  int pos;
  int base;
  double log_quality;
  bool reverse_complement;
  bool prior;

  VoteInfo(int pos, int base, double log_quality, bool reverse_complement, bool prior=false):pos(pos),base(base),log_quality(log_quality),reverse_complement(reverse_complement),prior(prior) {}
};


void graph_update(gl_types::iscope &scope,
    gl_types::icallback &scheduler);
#endif
