// Implementation for kinetic monte carlo engine

#include <cassert>

#include "kmc_engine.h"

// XXX hack to get hardcoded version working
#include "md_bead.h"
#include "md_kmc_bead.h"

void kmcEngine::Init(space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies, ParticleTracking *pTracking, long seed) {
  space_ = pSpace;
  species_ = pSpecies;
  ndim_ = space_->n_dim;
  nperiodic_ = space_->n_periodic;
  tracking_ = pTracking;
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
}

// Initialize the simples, nsimples, etc
void kmcEngine::InitMP() {
  nsimples_ = tracking_->GetNSimples();
  simples_ = tracking_->GetSimples();
}

// Step the kmc engine forward one
void kmcEngine::StepKMC() {
  PrepKMC();

  HardcodedBindUnbind();
  // Ask each species to do their own StepKMC
  /*for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->IsKMC()) {
      (*spec)->StepKMC();
    }
  }*/
}

// XXX CJE fix this later
void kmcEngine::HardcodedBindUnbind() {
  // Randomly choose order of binding and unbinding
  int g[2] = {0, 1};
  for (int i = 0; i < 2; ++i) {
    int j = gsl_rng_uniform_int(rng_.r, 2);
    int swapme = g[i];
    g[i] = g[j];
    g[j] = swapme;
  }

  for (int i = 0; i < 2; ++i) {
    switch (g[i]) {
      case 0:
        HardcodedBind();
        break;
      case 1:
        HardcodedUnbind();
        break;
    }
  }
}

// XXX CJE fix this later
void kmcEngine::HardcodedBind() {
  double eps_eff_ = 30910;
  double on_rate_ = 0.003916;
  // Loop over particles and get both if interacting
  for (int idx = 0; idx < nsimples_; ++idx) {
    auto part = (*simples_)[idx];
    if (!part->IsKMC()) continue;
    // Dynamic cast to MDKMCBead
    MDKMCBead *pkmcbead = dynamic_cast<MDKMCBead*>(part);
    double binding_affinity = eps_eff_ * on_rate_ * pkmcbead->GetDelta();
    auto nexp = pkmcbead->GetNExp();
    if (!pkmcbead->GetBound() && nexp > 0.0) {
      auto rng = pkmcbead->GetRNG();
      double roll = gsl_rng_uniform(rng->r);
      if (roll < nexp) {
        printf("[%d] Successful KMC move {nexp: %2.4f}, {roll: %2.4f}\n", pkmcbead->GetOID(), nexp, roll);
        pkmcbead->SetBound(true);
        pkmcbead->SetNExp(0.0);
        // Figure out where to attach
        double pos = 0.0;
        auto neighbors = tracking_->GetNeighbors();
        for (auto nldx = neighbors[idx].begin(); nldx != neighbors[idx].end(); ++nldx) {
          pos += binding_affinity * nldx->kmc_;
          if (pos > roll) {
            int jdx = nldx->idx_;
            auto part2 = (*simples_)[jdx];
            MDBead* pbead = dynamic_cast<MDBead*>(part2);
            printf("[%d,%d] Attaching to [%d,%d] {pos: %2.4f}\n", idx, pkmcbead->GetOID(), jdx, pbead->GetOID(), pos);
            pkmcbead->Attach(jdx);
            break;
          }
        } //which one to attach to
      } // successful kmc move
    } // are we already bound?

    // If we are bound, update position and velocity to the attached part
    if (pkmcbead->GetBound()) {
      pkmcbead->SetNExp(0.0);
      auto aidx = pkmcbead->GetAttach();
      auto part2 = (*simples_)[aidx];
      auto apos = part2->GetRigidPosition();
      auto apos_scaled = part2->GetRigidScaledPosition();
      auto vpos = part2->GetVelocity();

      auto mpos = pkmcbead->GetRigidPosition();
      printf("[%d,%d] attached [%d,%d], (%2.2f, %2.2f) -> setting (%2.2f, %2.2f)\n",
              idx, pkmcbead->GetOID(), aidx, part2->GetOID(), mpos[0], mpos[1],
              apos[0], apos[1]);
      pkmcbead->SetPrevPosition(mpos);
      pkmcbead->SetPosition(apos);
      pkmcbead->SetScaledPosition(apos_scaled);
      pkmcbead->SetVelocity(vpos);
    }
  }
}

void kmcEngine::HardcodedUnbind() {

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
  // Ask the species to update
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->IsKMC()) {
      (*spec)->PrepKMC();
    }
  }
}
