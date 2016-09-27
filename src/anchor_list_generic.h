#ifndef _SIMCORE_ANCHOR_LIST_GENERIC_H_
#define _SIMCORE_ANCHOR_LIST_GENERIC_H_

// Anchor points are for tethering potentials, when we have
// a particle attached to a specific positon, or other particle
// and have to update this.  Especially important for rigid objects
// like SPBs, centrosomes

#include <vector>
#include <unordered_map>

struct _anchor {
  int idx_base_;
  int idx_other_;
  double pos0_[3]; // lab position of first particle
  double pos_rel0_[3]; // relative position if on structure
  double pos1_[3]; // lab position of second particle
  double pos_rel1_[3];
};
typedef struct _anchor anchor_t;

typedef std::vector<anchor_t> al_list;
typedef std::unordered_map<int, std::vector<anchor_t>> al_set;

#endif
