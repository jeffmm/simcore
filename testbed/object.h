#ifndef _SIMCORE_OBJECT_H_
#define _SIMCORE_OBJECT_H_

#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include "rng.h"
#include "parameters.h"

class Object {
  private:
    static unsigned int _next_oid_;
    static unsigned int _n_dim_;
    static long _seed_;
  protected:
    RNG rng_;
    unsigned int n_dim_;
    unsigned int oid_;
    double position_[3],
           orientation_[3],
           length_,
           diameter_;
    system_parameters const * params_;
  public:
    Object();
    static void SetNDim(unsigned int n_dim) {_n_dim_ = n_dim;}
    static void SetSeed(long seed) {_seed_ = seed;}
    unsigned int GetOID();
    double const * const GetPosition();
    double const * const GetOrientation();
    double const GetLength();
    double const GetDiameter();
    void SetLength(double l);
    void SetDiameter(double d);
    void SetPosition(double * pos);
    void SetOrientation(double * u);
    virtual void Report();
};

class Bond; // Forward declaration
// Sites, ie graph vertices
class Site : public Object {
  protected:
    std::vector<Bond*> bonds_;
  public:
    Site() {}
    void AddBond(Bond * bond);
    void Report();
    Bond * GetBond(int i);
    Bond * GetOtherBond(unsigned int bond_oid);
};

// Bonds, ie graph edges
class Bond : public Object {
  protected:
    Site * sites_[2];
  public:
    Bond() {}
    Bond(Site * s1, Site * s2);
    void Init(Site * s1, Site * s2);
    void Report();
    Site * GetSite(int i);
    Bond * GetNeighborBond(int i);
};

class Mesh : public Object {
  protected:
    unsigned int n_sites_;
    unsigned int n_bonds_;
    std::vector<Site> sites_;
    std::vector<Bond> bonds_;
  public:
    Mesh() {n_sites_ = n_bonds_ = 0;}
    void InitSiteAt(double * pos, double d);
    void InitBondAt(double * pos, double * u, double l, double d);
    void InitRandomSite(double d);
    void AddRandomBondToSite(double l, int i_site);
    void AddRandomBondAnywhere(double l, double d);
    void AddRandomBondAtTip(double l, double d);
    void AddBondAtTip(double *u, double l, double d);
    void AddSite(Site s);
    void AddBond(Bond b);
    void ReportSites();
    void ReportBonds();
    void Report();
};

#endif
