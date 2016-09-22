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
    AddKMC(sid1, sid2, new_kmc);
  }
}

void kmcEngine::AddKMC(SID sid1, SID sid2, KMCBase *kmc) {
  sid_pair key1 = std::make_pair(sid1, sid2);
  if (kmc_map_.count(key1)) return;
  kmc_map_[key1] = kmc;
}

double kmcEngine::GetMaxRcut() {
  double max_rcut = 0.0;
  for (auto kmc = kmc_map_.begin(); kmc != kmc_map_.end(); ++kmc) {
    max_rcut = std::max(max_rcut, kmc->second->GetMaxRcut());
  }
  return max_rcut;
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
  // Move through the map and ask it to do the kmc prep
  for (auto kmc = kmc_map_.begin(); kmc != kmc_map_.end(); ++kmc) {
    if (debug_trace)
      printf("PrepKMC [%hhu,%hhu]\n", kmc->first.first, kmc->first.second);
    kmc->second->PrepKMC();
  }
}

// Step the kmc engine forward one
void kmcEngine::StepKMC() {
  // Run the KMC modules that we know aboout
  // Loop over species pairs - we know that some self interact, but
  // some interact between species, and it is keyed on who is #1
  // and who is #2
  for (auto kmc = kmc_map_.begin(); kmc != kmc_map_.end(); ++kmc) {
    if (debug_trace)
      printf("StepKMC [%hhu,%hhu]\n", kmc->first.first, kmc->first.second);
    kmc->second->StepKMC();
  }
}

// Ask each species to do their own update on the kmc step
void kmcEngine::UpdateKMC() {
  for (auto kmc = kmc_map_.begin(); kmc != kmc_map_.end(); ++kmc) {
    if (debug_trace)
      printf("UpdateKMC [%hhu,%hhu]\n", kmc->first.first, kmc->first.second);
    kmc->second->UpdateKMC();
  }
}

// Transfer the forces (if need by) by the KMC
void kmcEngine::TransferForces() {
  for (auto kmc = kmc_map_.begin(); kmc != kmc_map_.end(); ++kmc) {
    kmc->second->TransferForces();
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
  for (auto kmc = kmc_map_.begin(); kmc != kmc_map_.end(); ++kmc) {
    kmc->second->Dump();
  }
  #endif
}

void kmcEngine::PrepOutputs() {
  for (auto kmc = kmc_map_.begin(); kmc != kmc_map_.end(); ++kmc) {
    kmc->second->PrepOutputs(); 
  }
}

void kmcEngine::WriteOutputs(int istep) {
  for (auto kmc = kmc_map_.begin(); kmc != kmc_map_.end(); ++kmc) {
    kmc->second->WriteOutputs(istep); 
  }
}
