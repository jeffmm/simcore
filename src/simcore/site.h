#ifndef _SIMCORE_SITE_H_
#define _SIMCORE_SITE_H_

#include "object.h"

class Bond; // Forward declaration
enum directed_type {OUTGOING,INCOMING,NONE};
typedef std::pair<Bond*,directed_type> directed_bond;
typedef std::vector<directed_bond>::iterator db_iterator;

// Sites, ie graph vertices
class Site : public Object {
  protected:
    std::vector<directed_bond> bonds_;
    double tangent_[3]; // if one or two bonds, vector tangent to bonds at site
    double random_force_[3]; // random forces for filaments
    int n_bonds_;
  public:
    Site();
    void AddBond(Bond * bond, directed_type dir);
    void Report();
    void ReportBonds();
    Bond * GetBond(int i);
    Bond * GetOtherBond(int bond_oid);
    void CalcTangent();
    void SetRandomForce(double * f_rand);
    void AddRandomForce();
    double const * const GetRandomForce();
    double const * const GetTangent();
    directed_bond GetDirectedBond(int i);
    directed_bond GetOutgoingBond();
    directed_bond GetOtherDirectedBond(int bond_oid);
    void RemoveOutgoingBonds();
    void RemoveBond(int bond_oid);
    virtual bool HasNeighbor(int other_oid);
    void Draw(std::vector<graph_struct*> * graph_array);
};

#endif // _SIMCORE_SITE_H_
