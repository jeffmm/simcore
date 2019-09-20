#ifndef _SIMCORE_CROSSLINK_H_
#define _SIMCORE_CROSSLINK_H_

//#include "species.hpp"
#include "anchor.hpp"
#include "minimum_distance.hpp"
#include <kmc.hpp>
#include <mutex>

enum xstate { unbound, singly, doubly };

class Neighbors {
private:
  std::mutex mtx_;
  std::vector<Object *> nlist_;

public:
  Neighbors() {}
  ~Neighbors() { Clear(); }
  Neighbors(const Neighbors &that) { this->nlist_ = that.nlist_; }
  Neighbors &operator=(Neighbors const &that) {
    this->nlist_ = that.nlist_;
    return *this;
  }
  void AddNeighbor(Object *obj) {
    std::lock_guard<std::mutex> lk(mtx_);
    nlist_.push_back(obj);
  }
  const Object *const *GetNeighborsMem() { return &nlist_[0]; }
  void Clear() { nlist_.clear(); }
  int NNeighbors() { return nlist_.size(); }
  Object *GetNeighbor(int i_neighbor) {
    if (i_neighbor >= nlist_.size()) {
      error_exit("Invalid index received in class Neighbor");
    }
    return nlist_[i_neighbor];
  }
};

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
  Neighbors neighbors_;
  std::vector<int> kmc_filter_;
  void CalculateTetherForces();
  void AttemptCrosslink();
  void CalculateBinding();
  void CalcBinding();
  void SinglyKMC();
  void DoublyKMC();

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
  void ClearNeighbors();
};

#endif
