// All neighbor list implementations and information

#ifndef _SIMCORE_NEIGHBOR_LIST_H_
#define _SIMCORE_NEIGHBOR_LIST_H_

#include <vector>

struct _neighbor_kmc {
  int idx_;
  double kmc_ = 0.0;
};
typedef struct _neighbor_kmc neighbor_kmc_t;
typedef std::vector<neighbor_kmc_t> nl_kmc_list;

#endif
