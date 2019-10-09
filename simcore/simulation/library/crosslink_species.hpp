#ifndef _SIMCORE_CROSSLINK_SPECIES_H_
#define _SIMCORE_CROSSLINK_SPECIES_H_

#include "crosslink.hpp"
#include "species.hpp"
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

typedef std::vector<std::pair<std::vector<Crosslink>::iterator,
                              std::vector<Crosslink>::iterator>>
    xlink_chunk_vector;
typedef std::vector<Crosslink>::iterator xlink_iterator;

class CrosslinkSpecies : public Species<Crosslink, species_id::crosslink> {
private:
  bool *update_;
  std::string checkpoint_file_;
  double *obj_volume_;
  double xlink_concentration_;
  double k_on_;
  double k_off_;
  double k_on_d_;
  double k_off_d_;
  RNG rng_;
  LookupTable lut_;
  std::vector<Object *> *objs_;
  void CalculateBindingFree();
  void BindCrosslink();
  void UpdateBoundCrosslinks();
  void UpdateBoundCrosslinkForces();
  void UpdateBoundCrosslinkPositions();
  void ApplyCrosslinkTetherForces();
  Object *GetRandomObject();

public:
  void Init(system_parameters *params, species_base_parameters *sparams,
            space_struct *space);
  void InitInteractionEnvironment(std::vector<Object *> *objs, double *obj_vol,
                                  bool *update);
  void GetInteractors(std::vector<Object *> &ixors);
  void UpdatePositions();
  void CleanUp();
  void Draw(std::vector<graph_struct *> &graph_array);
  void BindCrosslinkObj(Object *obj);
  void AddNeighborToAnchor(Object *anchor, Object *neighbor);
  void AddMember();
  void GetAnchorInteractors(std::vector<Object *> &ixors);
  const double GetConcentration() const { return sparams_.concentration; }
  void ReadSpecs();

};

#endif
