#ifndef _SIMCORE_CROSSLINK_H_
#define _SIMCORE_CROSSLINK_H_

//#include "species.hpp"
#include "anchor.hpp"
#include "minimum_distance.hpp"
#include <kmc.hpp>
#include <kmc_choose.hpp>

enum xstate { unbound, singly, doubly };

/* Class that represents a two-headed crosslink that can create tethers between
 * two objects in the simulation. Governs binding and unbinding of crosslink
 * heads, and tether forces between bound heads. */
class Crosslink : public Object {
private:
  MinimumDistance *mindist_;
  draw_type draw_;
  bind_state state_;
  LookupTable *lut_;
  std::vector<int> kmc_filter_;
  double k_on_;
  double k_on_d_;
  double k_off_;
  double k_off_d_;
  double k_spring_;
  double k_align_;
  double f_spring_max_;
  double rest_length_;
  double rcapture_;
  double tether_force_;
  double fdep_factor_;
  double polar_affinity_;
  std::vector<Anchor> anchors_;
  void CalculateTetherForces();
  void CalculateBinding();
  void SinglyKMC();
  void DoublyKMC();
  void UpdateAnchorsToMesh();
  void UpdateAnchorPositions();
  void UpdateXlinkState();

public:
  Crosslink();
  void Init(crosslink_parameters *sparams, MinimumDistance *mindist,
            LookupTable *lut);
  void AttachObjRandom(Object *obj);
  void UpdateCrosslinkForces();
  void UpdateCrosslinkPositions();
  void GetAnchors(std::vector<Object *> &ixors);
  void GetInteractors(std::vector<Object *> &ixors);
  void Draw(std::vector<graph_struct *> *graph_array);
  void SetDoubly();
  void SetSingly();
  void SetUnbound();
  bool IsDoubly();
  bool IsUnbound();
  bool IsSingly();
  void UpdatePosition();
  void WriteSpec(std::fstream &ospec);
  void WriteCheckpoint(std::fstream &ocheck);
  void ReadSpec(std::fstream &ispec);
  void ReadCheckpoint(std::fstream &icheck);
  void ClearNeighbors();
  void ZeroForce();
  void ApplyTetherForces();
};

#endif
