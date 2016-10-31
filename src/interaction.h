#ifndef _SIMCORE_INTERACTION_H_
#define _SIMCORE_INTERACTION_H_

#include "anchor_list_generic.h"
#include "auxiliary.h"
#include "neighbor_list.h"

enum class ptype: int {
  external = 0,
  kmc = 1,
  internal = 2,
  tether = 3,
  boundary = 4,
  NONE = 5
};

inline ptype StringToPtype(const std::string& s) {
  if (s == "external")
    return ptype::external;
  else if (s == "kmc")
    return ptype::kmc;
  else if (s == "internal")
    return ptype::kmc;
  else if (s == "tether")
    return ptype::tether;
  else if (s == "boundary")
    return ptype::boundary;
  else {
    std::cout << "Wrong kind of potential interaction " << s << std::endl;
    exit(1);
    return ptype::NONE;
  }
}

inline std::string PtypeToString(const ptype& p) {
  std::string retval = "";
  if (p == ptype::external)
    retval = "external";
  else if (p == ptype::kmc)
    retval = "kmc";
  else if (p == ptype::internal)
    retval = "internal";
  else if (p == ptype::tether)
    retval = "tether";
  else if (p == ptype::boundary)
    retval = "boundary";
  return retval;
}

struct _interaction {
  int idx_;
  int jdx_;
  PotentialBase *pot_;
  ptype type_;
  int kmc_track_module_;

  // KMC stuff
  SID kmc_target_ = SID::none;
  neighbor_kmc_t *neighbor_ = nullptr;

  // Anchor
  anchor_t *anchor_ = nullptr;
};
typedef struct _interaction interaction_t;

#endif // _SIMCORE_INTERACTION_H_
