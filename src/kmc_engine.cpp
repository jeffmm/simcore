// Implementation for kinetic monte carlo engine

#include <cassert>

#include "kmc_md_mdkmc_bindunbind.h"
#include "kmc_engine.h"

// XXX hack to get hardcoded version working
#include "md_bead.h"
#include "md_kmc_bead.h"

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

// Register the available KMC modules
void kmcEngine::RegisterKMC() {
  REGISTER_KMC(MdMdkmcBindUnbind);
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

  nkmcs_ = node["kmc"].size();

  for (int ikmc = 0; ikmc < nkmcs_; ++ikmc) {
    std::string kmcname = node["kmc"][ikmc]["type"].as<std::string>();
    std::string sid1s   = node["kmc"][ikmc]["sid1"].as<std::string>();
    std::string sid2s   = node["kmc"][ikmc]["sid2"].as<std::string>();
    SID sid1 = StringToSID(sid1s);
    SID sid2 = StringToSID(sid2s);
    KMCBase *new_kmc = (KMCBase*) kmc_factory_.construct(kmcname);
    new_kmc->Init(space_, tracking_, ikmc, node, gsl_rng_get(rng_.r));
    AddKMC(sid1, sid2, new_kmc);
  }
}

void kmcEngine::AddKMC(SID sid1, SID sid2, KMCBase *kmc) {
  sid_pair key1 = std::make_pair(sid1, sid2);
  if (kmc_map_.count(key1)) return;
  kmc_map_[key1] = kmc;
}

// Run the overall kmc routine one step
void kmcEngine::RunKMC() {
  // Check to see if an update was triggered
  if (tracking_->TriggerUpdate()) {
    nsimples_ = tracking_->GetNSimples();
    simples_ = tracking_->GetSimples();
  }
  // Prepare the kmc probabilities, etc
  PrepKMC();

  // Run through our modules
  StepKMC();

  // Ask all the species to do their own updating
  UpdateKMC();
}

// Prepare and update probabilities of the kmc engine
void kmcEngine::PrepKMC() {
  // Have to go through the neighbor list and ask what things are doing...
  auto neighbors = tracking_->GetNeighbors();
  for (int idx = 0; idx < nsimples_; ++idx) {
    auto part = (*simples_)[idx];
    if (!part->IsKMC()) continue;
    part->PrepKMC(&neighbors[idx]);
  }
  // Ask the species to prepare
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->IsKMC()) {
      (*spec)->PrepKMC();
    }
  }
}

// Step the kmc engine forward one
void kmcEngine::StepKMC() {
  // Run the KMC modules that we know aboout
  // Loop over species pairs - we know that some self interact, but
  // some interact between species, and it is keyed on who is #1
  // and who is #2
  for (auto spec1 = species_->begin(); spec1 != species_->end(); ++spec1) {
    for (auto spec2 = species_->begin(); spec2 != species_->end(); ++spec2) {
      SID sid1 = (*spec1)->GetSID();
      SID sid2 = (*spec2)->GetSID();

      KMCBase *kmc = kmc_map_[std::make_pair(sid1, sid2)];
      if (kmc == nullptr) continue;
      kmc->RunKMC(*spec1, *spec2);
    }
  }
}

// Ask each species to do their own update on the kmc step
void kmcEngine::UpdateKMC() {
  // Ask each species to do their own update
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->IsKMC()) {
      (*spec)->StepKMC();
    }
  }
}

void kmcEngine::Print() {
  printf("********\n");
  printf("kmcEngine -> print\n");
  for (auto kmc = kmc_map_.begin(); kmc != kmc_map_.end(); ++kmc) {
    printf("{%d,%d} : ", (int)kmc->first.first, (int)kmc->first.second);
    kmc->second->Print();
  }
}

void kmcEngine::Dump() {
  #ifdef DEBUG
  printf("--------\n");
  printf("kmcEngine -> dump\n");
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->IsKMC()) {
      (*spec)->DumpKMC();
    }
  }
  #endif
}
