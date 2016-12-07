// Implementation for kinetic monte carlo engine

#include <cassert>

#include "dynamic_instability.h"
#include "xlink_kmc.h"
#include "kmc_engine.h"

#define REGISTER_KMC(n) kmc_factory_.register_class<n>(#n);

void kmcEngine::Init(space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, ParticleTracking *pTracking, long seed, char *pFname) {
  space_ = pSpace;
  species_ = pSpecies;
  ndim_ = space_->n_dim;
  nperiodic_ = space_->n_periodic;
  tracking_ = pTracking;
  rng_.init(seed);
  fname_ = pFname;
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
}

void kmcEngine::InitPotentials(PotentialManager *pPotentials) {
  potentials_ = pPotentials;
}

// Register the available KMC modules
void kmcEngine::RegisterKMC() {
  REGISTER_KMC(XlinkKMC);
  REGISTER_KMC(DynamicInstabilityKMC);
}

// Initialize the simples, nsimples, etc
void kmcEngine::InitMP() {
  nsimples_ = tracking_->GetNSimples();
  simples_ = tracking_->GetSimples();

  ParseKMC();
}

void kmcEngine::ParseKMC() {
  // Initialize based on fname
  std::cout << "********\n";
  std::cout << "KMC Load ->\n";
  std::cout << "  file: " << fname_ << std::endl;
  YAML::Node node = YAML::LoadFile(fname_);

  nkmcs_ = (int)node["kmc"].size();

  kmc_modules_.clear();
  kmc_modules_.resize(nkmcs_);

  for (int ikmc = 0; ikmc < nkmcs_; ++ikmc) {
    std::string kmcname = node["kmc"][ikmc]["type"].as<std::string>();
    std::string sid1s   = node["kmc"][ikmc]["sid1"].as<std::string>();
    std::string sid2s   = node["kmc"][ikmc]["sid2"].as<std::string>();
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
    KMCBase *new_kmc = (KMCBase*) kmc_factory_.construct(kmcname);
    new_kmc->Init(space_, tracking_, potentials_, spec1, spec2, ikmc, node, gsl_rng_get(rng_.r));
    kmc_modules_[ikmc] = new_kmc;
  }
}

double kmcEngine::GetMaxRcut() {
  double max_rcut = 0.0;
  for (auto kmc = kmc_modules_.begin(); kmc != kmc_modules_.end(); ++kmc) {
    max_rcut = std::max(max_rcut, (*kmc)->GetMaxRcut());
  }
  return max_rcut;
}

// Run the modules IN ORDER!!!!
void kmcEngine::RunKMC() {
  for (auto kmc = kmc_modules_.begin(); kmc != kmc_modules_.end(); ++kmc) {
    (*kmc)->PrepKMC();
    (*kmc)->StepKMC();
    (*kmc)->UpdateKMC();
    (*kmc)->TransferForces();
  }
}

void kmcEngine::Print() {
  std::cout << "********\n";
  std::cout << "kmcEngine -> print\n";
  for (auto kmc = kmc_modules_.begin(); kmc != kmc_modules_.end(); ++kmc) {
    auto msids = (*kmc)->GetSIDs();
    std::cout << "{" << (int)msids.first << "," << (int)msids.second << "} : ";
    (*kmc)->Print();
  }
}

void kmcEngine::Dump() {
  #ifdef DEBUG
  std::cout << "--------\n";
  std::cout << "kmcEngine -> dump\n";
  for (auto kmc = kmc_modules_.begin(); kmc != kmc_modules_.end(); ++kmc) {
    (*kmc)->Dump();
  }
  #endif
}

void kmcEngine::PrepOutputs() {
  for (auto kmc = kmc_modules_.begin(); kmc != kmc_modules_.end(); ++kmc) {
    (*kmc)->PrepOutputs(); 
  }
}

void kmcEngine::WriteOutputs(int istep) {
  for (auto kmc = kmc_modules_.begin(); kmc != kmc_modules_.end(); ++kmc) {
    (*kmc)->WriteOutputs(istep); 
  }
}
