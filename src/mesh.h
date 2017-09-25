#ifndef _SIMCORE_MESH_H_
#define _SIMCORE_MESH_H_

#include "bond.h"

typedef std::vector<Bond>::iterator bond_iterator;
typedef std::vector<Site>::iterator site_iterator;

class Mesh : public Object {
  private:
    static int next_mesh_id_;
  protected:
    int n_sites_,
        n_bonds_,
        n_bonds_max_;
    std::vector<Site> sites_;
    std::vector<Bond> bonds_;
    double bond_length_;
    Bond * GetRandomBond();
  public:
    Mesh();
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
    void RemoveBondFromTip();
    void ReportSites();
    void ReportBonds();
    void Report();
    void SubReport();
    void UpdateBondPositions();
    void UpdatePrevPositions();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    void Reserve(int n_bonds);
    void Clear();
    void DoubleGranularityLinear();
    void HalfGranularityLinear();
    Site * GetSite(int i);
    Bond * GetBond(int i);
    virtual void ZeroForce();
    std::vector<Object*> GetInteractors();
    virtual int GetCount();
    virtual void ReadPosit(std::fstream &ip);
    virtual void WritePosit(std::fstream &op);
    virtual void ReadSpec(std::fstream &ip);
    virtual void WriteSpec(std::fstream &op);
    virtual void ReadCheckpoint(std::fstream &ip);
    virtual void WriteCheckpoint(std::fstream &op);
    virtual void ScalePosition();
    virtual void UpdateDrTot();
    virtual double const GetDrTot();
    virtual void ZeroDrTot();
};

#endif // _SIMCORE_MESH_H_
