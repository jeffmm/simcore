#include "crosslink_manager.hpp"

void CrosslinkManager::Init(system_parameters *params, space_struct *space,
                            MinimumDistance *mindist,
                            std::vector<Object *> *objs) {
  objs_ = objs;
  update_ = false;
  obj_volume_ = 0;
  params_ = params;
  space_ = space;
  mindist_ = mindist;
}

void CrosslinkManager::InitSpecies(sid_label &slab, ParamsParser &parser) {
  xlink_species_.push_back(new CrosslinkSpecies());
  species_base_parameters *sparams = parser.GetNewSpeciesParameters(slab);
  xlink_species_.back()->Init(params_, sparams, space_);
  delete sparams;
  xlink_species_.back()->InitInteractionEnvironment(mindist, objs, obj_volume_,
                                                    update_);
}

/* Keep track of volume of objects in the system. Affects the
 * probability of a free crosslink binding to an object. */
void CrosslinkManager::UpdateObjsVolume() {
  obj_volume_ = 0;
  for (auto it = objs_->begin(); it != objs_->end(); ++it) {
    obj_volume_ += (*it)->GetVolume();
  }
}

/* Whether to reinsert anchors into the interactors list */
bool CrosslinkManager::CheckUpdate() {
  if (update_) {
    update_ = false;
    return true;
  }
  return false;
}

/* Return singly-bound anchors, for finding neighbors to bind to */
void CrosslinkManager::GetInteractors(std::vector<Object *> &ixors) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->GetInteractors(ixors);
  }
}

/* Returns all anchors, not just singly-bound anchors. Used for reassigning
   bound anchors to bonds upon a checkpoint reload */
void CrosslinkManager::GetAnchorInteractors(std::vector<Object *> &ixors) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->GetAnchorInteractors(ixors);
  }
}

void CrosslinkManager::UpdateCrosslinks() {
  update_ = false;
  UpdateObjsVolume();
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->UpdateCrosslinks();
  }
}

void CrosslinkManager::Clear() {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->CleanUp();
    delete (*it);
  }
  xlink_species_.clear();
}

void CrosslinkManager::Draw(std::vector<graph_struct *> *graph_array) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->Draw(graph_array);
  }
}

void CrosslinkManager::AddNeighborToAnchor(Object *anchor, Object *neighbor) {
  if (anchor->GetSID() != +species_id::crosslink) {
    Logger::Error(
        "AddNeighborToAnchor expected crosslink object, got generic object.");
  }
  Anchor *a = dynamic_cast<Anchor *>(anchor);
  a->AddNeighbor(neighbor);
}

void CrosslinkManager::InitOutputs(bool reading_inputs, bool reduce_flag,
                                   bool with_reloads) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->InitOutputs(reading_inputs, reduce_flag, with_reloads);
  }
}

void CrosslinkManager::WriteOutputs() {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->WriteOutputs();
  }
}

void CrosslinkManager::ReadInputs() {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->ReadInputs();
  }
}

void CrosslinkManager::InitCrosslinkSpecies() {}
