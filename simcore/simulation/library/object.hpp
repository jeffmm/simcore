#ifndef _SIMCORE_OBJECT_H_
#define _SIMCORE_OBJECT_H_

#include "auxiliary.hpp"
#include "interaction.hpp"
#include "rng.hpp"
#include <mutex>

class Object {
private:
  int oid_;
  int mesh_id_;
  static int _next_oid_;
  static std::mutex _obj_mtx_;
  void InitOID();

protected:
  static system_parameters *params_;
  static space_struct *space_;
  static int n_dim_;
  static double delta_;
  species_id sid_;
  obj_type type_;
  graph_struct g_;
  RNG rng_;
  draw_type draw_;
  int n_contact_;
  int in_flock_;           // 0 if not in flock, 1 if interior, 2 if exterior
  int flock_change_state_; // 0 if same as previous step, 1 if joined flock, 2
                           // if left flock
  double position_[3];
  double prev_position_[3];
  double prev_orientation_[3];
  double scaled_position_[3];
  double orientation_[3];
  double force_[3];
  double torque_[3];
  double dr_zero_[3];
  double color_;
  double diameter_;
  double length_;
  double p_energy_;
  double dr_tot_;
  double polar_order_;
  double contact_number_;
  bool interacting_;
  bool is_mesh_;
  bool has_overlap_;
  bool interactor_update_;

  std::vector<Object *> interactors_;
  std::vector<Interaction *> ixs_;
  void UpdateKMC();

public:
  Object();
  // kmcx parameter
  int gid;
  double length;
  double radius;
  double pos[3];
  double direction[3];

  // Static functions
  static void SetParams(system_parameters *params);
  static void SetSpace(space_struct *space);
  static void SetNDim(int n_dim);
  static void SetDelta(double delta);
  // Trivial Set/Get functions
  void SetSID(species_id sid);
  void SetType(obj_type type);
  void SetPosition(double const *const pos);
  void SetScaledPosition(double const *const spos);
  void SetOrientation(double const *const u);
  void SetPrevPosition(double const *const ppos);
  void SetPrevOrientation(double const *const pu);
  void SetDiameter(double new_diameter);
  void SetLength(double new_length);
  void AddForce(double const *const f);
  void SubForce(double const *const f);
  void SetForce(double const *const f);
  void AddTorque(double const *const t);
  void SubTorque(double const *const t);
  void SetTorque(double const *const t);
  void AddPotential(double const p);
  void AddPolarOrder(double const po);
  void AddContactNumber(double const cn);
  void SetInteractor(bool ix);
  void ToggleIsMesh();
  void CalcPolarOrder();
  void ZeroPolarOrder();
  species_id const GetSID();
  obj_type const GetType();
  int const GetOID() const;
  int const GetMeshID() const;
  double const *const GetPosition();
  double const *const GetPrevPosition();
  double const *const GetScaledPosition();
  double const *const GetOrientation();
  virtual void GetAvgPosition(double *ap);
  virtual void GetAvgOrientation(double *au);
  virtual void SetAvgPosition();
  double const GetDiameter();
  double const GetLength();
  double const *const GetForce();
  double const *const GetTorque();
  double const GetPotentialEnergy();
  double const GetPolarOrder();
  double const GetContactNumber();
  bool const IsInteractor();
  bool const IsMesh();
  bool const CheckInteractorUpdate();
  void HasOverlap(bool overlap);
  void SetFlockType(int in_flock);
  void SetFlockChangeState(int fcs);
  int GetFlockType();
  int GetFlockChangeState();
  void SetMeshID(int mid);

  // Virtual functions
  virtual void Init(species_base_parameters *sparams) {}
  virtual void InsertRandom();
  virtual void InsertRandomOriented(double *u);
  virtual void InsertAt(double *pos, double *u);
  virtual void ZeroForce();
  virtual void UpdatePeriodic();
  virtual void UpdatePosition() {}
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  virtual void SetColor(double const c, draw_type dtype);
  virtual void ScalePosition();
  virtual int GetCount();
  virtual void GetInteractors(std::vector<Object *> &ix);
  virtual double const *const GetInteractorPosition();
  virtual double const *const GetInteractorPrevPosition();
  virtual double const *const GetInteractorScaledPosition();
  virtual double const *const GetInteractorOrientation();
  virtual double const GetInteractorDiameter();
  virtual double const GetInteractorLength();
  virtual double const GetVolume();
  virtual void UpdateDrTot();
  virtual double const GetDrTot();
  virtual void ZeroDrTot();
  virtual bool HasNeighbor(int other_id);
  virtual void GiveInteraction(Interaction *ix);
  virtual std::vector<Interaction *> *GetInteractions();
  virtual void ClearInteractions();
  virtual void Cleanup();
  // virtual void BindAnchor(anchor *ix);
  // virtual void UnbindAnchor();

  // I/O functions
  virtual void Report();
  virtual void WritePosit(std::fstream &oposit);
  virtual void ReadPosit(std::fstream &iposit);
  virtual void WriteSpec(std::fstream &ospec);
  virtual void ReadSpec(std::fstream &ispec);
  virtual void ReadPositFromSpec(std::fstream &ispec);
  virtual void WriteCheckpoint(std::fstream &ocheck);
  virtual void ReadCheckpoint(std::fstream &icheck);
  void SetRNGState(const std::string &filename);
};

// void MinimumDistance(Object* o1, Object* o2, Interaction *ix, space_struct
// *space); void BoundaryConditions(Object * o1, space_struct *space);

#endif // _SIMCORE_OBJECT_H_
