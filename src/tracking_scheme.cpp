// Implementation of tracking scheme base

#include "tracking_scheme.h"

// Initialize
void TrackingScheme::Init(space_struct *pSpace,
                          PotentialBase *pPotentialBase,
                          std::vector<interaction_t> *pInteractions,
                          YAML::Node *pNode) {

  space_ = pSpace;
  pbase_ = pPotentialBase;
  ndim_ = space_->n_dim;
  nperiodic_ = space_->n_periodic;
  interactions_ = pInteractions;

  #ifdef ENABLE_OPENMP
  #pragma omp parallel
  {
    if (0 == omp_get_thread_num()) {
      nthreads_ = omp_get_num_threads();
    }
  }
  #else
  nthreads_ = 1;
  #endif

  YAML::Node node = *pNode;

  // SIDs for this interaction
  std::string sid0s = node["sid1"].as<std::string>();
  std::string sid1s = node["sid2"].as<std::string>();
  sid0_ = StringToSID(sid0s);
  sid1_ = StringToSID(sid1s);

  // Type of interaction
  std::string types = node["type"].as<std::string>();
  type_ = StringToPtype(types);
}

void TrackingScheme::Print() {
  std::cout << name_ << " Tracking Scheme" << std::endl;
  std::cout << "   potential: ->\n";
  std::cout << "   {" << SIDToString(sid0_) << ", " << SIDToString(sid1_) << "} : ";
  pbase_->Print();
  std::cout << "   type: " << PtypeToString(type_) << std::endl;
}
