#ifndef _SIMCORE_SITE_H_
#define _SIMCORE_SITE_H_

#include "object.h"

class Bond; // Forward declaration
enum directed_type {OUTGOING,INCOMING,NONE};
typedef std::pair<Bond*,directed_type> directed_bond;

// Sites, ie graph vertices
class Site : public Object {
  protected:
    std::vector<directed_bond> bonds_;
    int n_bonds_;
  public:
    Site() : Object() {n_bonds_ = 0;}
    void AddBond(Bond * bond, directed_type dir);
    void Report();
    void ReportBonds();
    Bond * GetBond(int i);
    Bond * GetOtherBond(int bond_oid);
    directed_bond GetOtherDirectedBond(int bond_oid);

};

#endif
