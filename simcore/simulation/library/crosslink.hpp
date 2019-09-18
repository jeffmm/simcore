#ifndef _SIMCORE_CROSSLINK_H_
#define _SIMCORE_CROSSLINK_H_

//#include "species.hpp"
#include "anchor.hpp"
#include "minimum_distance.hpp"
#include <kmc.hpp>
//#include <mutex>

enum xstate { unbound, singly, doubly };

class Crosslink : public Object {
private:
  Interaction ix;
  MinimumDistance *mindist_;
  draw_type draw_;
  bind_state state_;
  LookupTable *lut_;
  double k_on_;
  double k_on_sd_;
  double k_off_;
  double k_spring_;
  double k_align_;
  double f_spring_max_, rest_length_;
  double rcapture_;
  double tether_force_;
  double fdep_factor_;
  double polar_affinity_;
  std::vector<Anchor> anchors_;
  std::vector<Object *> nlist_;
  std::vector<int> kmc_filter_;
  void CalculateTetherForces();
  void AttemptCrosslink();
  void CalculateBinding();
  void CalcBinding();
  void SinglyKMC();
  void DoublyKMC();
  /* TODO, get rid of racy neighborlist additions */
  // std::mutex xlink_mtx_;

public:
  Crosslink();
  void Init(MinimumDistance *mindist, LookupTable *lut);
  void UnbindAnchor(bool second = false);
  void AttachObjRandom(Object *obj);
  void UpdateCrosslink();
  Anchor *GetBoundPtr();
  void GetAnchors(std::vector<Object *> *ixors);
  void Draw(std::vector<graph_struct *> *graph_array);
  void SetDoubly();
  void SetSingly();
  void SetUnbound();
  bool IsDoubly();
  bool IsUnbound();
  bool IsSingly();
  void UpdatePosition();
  void AddNeighbor(Object *neighbor);
  void WriteSpec(std::fstream &ospec);
  void WriteCheckpoint(std::fstream &ocheck);
  void ReadSpec(std::fstream &ispec);
  void ReadCheckpoint(std::fstream &icheck);
};

#endif
