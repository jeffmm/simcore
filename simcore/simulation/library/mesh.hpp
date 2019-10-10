#ifndef _SIMCORE_MESH_H_
#define _SIMCORE_MESH_H_

#include "bond.hpp"

typedef std::vector<Bond>::iterator bond_iterator;
typedef std::vector<Site>::iterator site_iterator;

class Mesh : public Object {
private:
  static int _next_mesh_id_;
  static std::mutex _mesh_mtx_;
  void InitMeshID();

protected:
  bool anchored_;
  bool midstep_;
  bool posits_only_;
  int n_sites_;
  int n_bonds_;
  int n_bonds_max_;
  std::vector<Site> sites_;
  std::vector<Bond> bonds_;
  double bond_length_;
  double true_length_ = -1;
  Bond *GetRandomBond();
  void UpdateInteractors();
  void UpdateSiteOrientations();
  void RelocateMesh(double *pos, double *u);
  void AddRandomBondToSite(double l, int i_site);
  void AddBondToSite(double *u, double l, int i_site);
  void AddSite(Site s);
  void AddBond(Site *s1, Site *s2);

public:
  Mesh();
  void InitSiteAt(double *pos, double d);
  void InitBondAt(double *pos, double *u, double l, double d);
  void InitRandomSite(double d);
  void AddRandomBondAnywhere(double l, double d = 1);
  void AddRandomBondToTip(double l);
  void AddBondToTip(double *u, double l);
  void SetBondLength(double l);
  void RemoveBondFromTip();
  void ReportSites();
  void ReportBonds();
  void Report();
  void SubReport();
  void UpdateBondPositions();
  void UpdatePrevPositions();
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  void Reserve(int n_bonds);
  void Clear();
  void DoubleGranularityLinear();
  void HalfGranularityLinear();
  int GetNBonds() { return n_bonds_; }
  Bond *GetBondAtLambda(double lambda);
  Site *GetSite(int i);
  Bond *GetBond(int i);
  virtual void ZeroForce();
  virtual void GetInteractors(std::vector<Object *> &ix);
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
  virtual void SetPosition(double const *const pos);
  virtual std::vector<Interaction *> *GetInteractions();
  virtual void ClearInteractions();
  virtual void GetAvgPosition(double *ap);
  virtual void GetAvgOrientation(double *au);
  virtual void SetAvgPosition();
  virtual void GetContactNumbers(std::vector<double> *cn);
  virtual void GetPolarOrders(std::vector<double> *po);
  std::pair<double, double> GetAvgOrientationCorrelation();
  virtual void ZeroOrientationCorrelations();
  virtual const double GetBondLength() const;
  virtual const bool CheckInteractorUpdate();
  virtual const double GetLambdaAtBond(int bond_oid);
  virtual const double GetTrueLength() const;
};

#endif // _SIMCORE_MESH_H_
