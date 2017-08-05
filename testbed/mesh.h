#include "bond.h"

class Mesh : public Object {
  protected:
    int n_sites_,
        n_bonds_,
        n_bonds_max_;
    std::vector<Site> sites_;
    std::vector<Bond> bonds_;
    double bond_length_;
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
    void ReportSites();
    void ReportBonds();
    void Report();
    void SubReport();
    void Draw(std::vector<graph_struct*> * graph_array);
    void Reserve(int n_bonds);
    Site * GetSite(int i);
    Bond * GetBond(int i);
};

