// All neighbor list implementations and information

#ifndef _SIMCORE_NEIGHBOR_LIST_H_
#define _SIMCORE_NEIGHBOR_LIST_H_

struct _neighbor_kmc {
  int idx_;
  double kmc_ = 0.0;
};
typedef struct _neighbor_kmc neighbor_kmc_t;

struct _neighbor_v2 {
  int oidx_;
  int ojdx_;
};
typedef struct _neighbor_v2 neighbor_v2_t;

#endif
