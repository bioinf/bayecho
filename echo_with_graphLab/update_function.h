#ifndef UPDATE_FUNCTION_
#define UPDATE_FUNCTION_
#include <graphlab.hpp>
#include <tr1/memory>
#include <graphlab/macros_def.hpp>
#include "node.h"
#include "edge.h"
typedef graphlab::graph<tr1::shared_ptr<MyNode>, tr1::shared_ptr<Edge> > graph_type;
typedef graphlab::types<graph_type> gl_types;

extern gl_types::glshared<double***> CONF_MAT;

void add_to_confMat(graphlab::any& current_confMat, const graphlab::any& added_value);

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
