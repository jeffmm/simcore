// Implementation of KMC engine v2

#include "kmc_engine_v2.h"

#include "xlink_kmc_v2.h"

#define REGISTER_KMC(n) kmc_factory_.register_class<n>(#n);

// Initialize with everything
void kmcEngineV2::Init(system_parameters *pParams,
                       space_struct *pSpace,
                       ParticleEngine *pTrackEngine,
                       std::vector<interaction_t> *pInteractions,
                       long seed) {
  params_ = pParams;
  space_ = pSpace;
  ptrack_ = pTrackEngine;
  interactions_ = pInteractions;

  ndim_ = space_->n_dim;
  nperiodic_ = space_->n_periodic;
  rng_.init(seed);

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

  RegisterKMC();

  AttachParticleEngine();
}

// Register the KMC modules
void kmcEngineV2::RegisterKMC() {
  REGISTER_KMC(XlinkKMCV2);
}

// Attach to the particle engine
void kmcEngineV2::AttachParticleEngine() {
  simples_ = ptrack_->GetSimples();
  species_ = ptrack_->GetSpecies();
}

// Initialize the MP stuff (really just get particle numbers)
void kmcEngineV2::InitMP() {
  nsimples_ = (int)simples_->size();
  nspecies_ = (int)species_->size();

  ParseKMC();
}

// Parse the KMC file
void kmcEngineV2::ParseKMC() {
  // Initialize based on file given to from params
  std::cout << "********\n";
  std::cout << "KMC Load ->\n";
  std::cout << "   file: " << params_->kmcfile << std::endl;

  YAML::Node node = YAML::LoadFile(params_->kmcfile);

  nkmcs_ = (int)node["kmc"].size();

  kmc_modules_.clear();
  kmc_modules_.resize(nkmcs_);

  for (int ikmc = 0; ikmc < nkmcs_; ++ikmc) {
    YAML::Node subnode = node["kmc"][ikmc];
    std::string kmcname = subnode["type"].as<std::string>();
    std::string sid1s   = subnode["sid1"].as<std::string>();
    std::string sid2s   = subnode["sid2"].as<std::string>();
    SID sid1 = StringToSID(sid1s);
    SID sid2 = StringToSID(sid2s);
    // Find the 2 species
    SpeciesBase *spec1 = nullptr, *spec2 = nullptr;
    for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
      if ((*spec)->GetSID() == sid1)
        spec1 = (*spec);
      if ((*spec)->GetSID() == sid2)
        spec2 = (*spec);
    }
    KMCBaseV2 *new_kmc = (KMCBaseV2*) kmc_factory_.construct(kmcname);
    new_kmc->Init(space_, ptrack_, spec1, spec2, &subnode, gsl_rng_get(rng_.r));
    kmc_modules_[ikmc] = new_kmc;
  }
}

// Run the modules IN ORDER!!!!
void kmcEngineV2::StepKMC() {
  // Run this in order
  for (auto kmc = kmc_modules_.begin(); kmc != kmc_modules_.end(); ++kmc) {
    (*kmc)->GenerateKMCNeighborList();
    (*kmc)->PrepKMC();
    (*kmc)->StepKMC();
    (*kmc)->UpdateKMC();
    (*kmc)->TransferForces();
  }
}

// Print stuff
void kmcEngineV2::Print() {
  std::cout << "********\n";
  std::cout << "kmcEngine -> print\n";
  for (auto kmc = kmc_modules_.begin(); kmc != kmc_modules_.end(); ++kmc) {
    auto msids = (*kmc)->GetSIDs();
    std::cout << "{" << (int)msids.first << "," << (int)msids.second << "} : ";
    (*kmc)->Print();
  }
}

// Dump stuff
void kmcEngineV2::Dump() {
  #ifdef DEBUG
  std::cout << "--------\n";
  std::cout << "kmcEngine -> dump\n";
  for (auto kmc = kmc_modules_.begin(); kmc != kmc_modules_.end(); ++kmc) {
    (*kmc)->Dump();
  }
  #endif
}
