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
    static int _next_oid_;
    static int _n_dim_;
    static long _seed_;
    static system_parameters const * _params_;
  protected:
    system_parameters const * params_;
    RNG rng_;
    int n_dim_;
    int oid_;
    double position_[3],
           scaled_position_[3],
           prev_position_[3],
           orientation_[3],
           force_[3],
           length_,
           diameter_;
    graph_struct g_struct;
  public:
    Object();
    static void SetNDim(int n_dim) {_n_dim_ = n_dim;}
    static void SetSeed(long seed) {_seed_ = seed;}
    static void SetParameters(system_parameters * params) {_params_=params;}
    int GetOID();
    double const * const GetPosition();
    double const * const GetOrientation();
    double const GetLength();
    double const GetDiameter();
    void SetLength(double l);
    void SetDiameter(double d);
    void SetPosition(double * pos);
    void SetOrientation(double * u);
    virtual void UpdatePosition();
    virtual void UpdatePrevPosition();
    virtual void Report();
    // Main draw function, return struct of graphics info
    virtual void Draw(std::vector<graph_struct*> * graph_array);
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
    Bond * GetOtherBond(int bond_oid);
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
    int n_sites_;
    int n_bonds_;
    std::vector<Site> sites_;
    std::vector<Bond> bonds_;
    double bond_length_=1;
  public:
    Mesh() {n_sites_ = n_bonds_ = 0;}
    void InitSiteAt(double * pos, double d);
    void InitBondAt(double * pos, double * u, double l, double d);
    void InitRandomSite(double d);
    void AddRandomBondToSite(double l, int i_site);
    void AddRandomBondAnywhere(double l, double d=1);
    void AddRandomBondToTip(double l);
    void AddBondToTip(double *u, double l);
    void AddBondToSite(double *u, double l, int i_site);
    void AddSite(Site s);
    void AddBond(Bond b);
    void SetBondLength(double l);
    void ReportSites();
    void ReportBonds();
    void Report();
    void Draw(std::vector<graph_struct*> * graph_array);
};

class Motor : public Site {

  protected:
    bool bound_;
    double k_on_,
           k_off_;
  public:
    Motor();
    bool UpdatePriors();
    void UpdatePosition();
    void Diffuse();
    void DiffuseAlongBond();
}

#endif
