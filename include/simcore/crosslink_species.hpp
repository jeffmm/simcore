#ifndef _SIMCORE_CROSSLINK_SPECIES_H_
#define _SIMCORE_CROSSLINK_SPECIES_H_

#include "crosslink.hpp"
#include "species.hpp"

typedef std::vector<std::pair<std::vector<Crosslink>::iterator,
                              std::vector<Crosslink>::iterator>>
    xlink_chunk_vector;
typedef std::vector<Crosslink>::iterator xlink_iterator;

class CrosslinkSpecies : public Species<Crosslink, species_id::crosslink> {
 private:
  bool *update_;
  std::string checkpoint_file_;
  double *obj_volume_;  // Total length of all the objects in the system
  double xlink_concentration_;
  double bind_site_density_;
  int begin_with_crosslinks_;
  bool infinite_reservoir_flag_;
  double k_on_;
  bool static_flag_;
  LookupTable lut_;
  std::vector<Object *> *objs_;

  LUTFiller *MakeLUTFiller();
  void CalculateBindingFree();
  void BindCrosslink();
  void UpdateBoundCrosslinks();
  void UpdateBoundCrosslinkForces();
  void UpdateBoundCrosslinkPositions();
  void ApplyCrosslinkTetherForces();
  Object *GetRandomObject();

 public:
  CrosslinkSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
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
  void ReadSpecs();
  void InsertCrosslinks();
  void BeginWithCrosslinks();
  void BeginWithCrosslinkers();
  const int GetDoublyBoundCrosslinkNumber() const;
  const double GetConcentration() const;
  const double GetRCutoff() const;
};

#endif
