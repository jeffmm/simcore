// Implementation for simple binding and unbinding

#include "kmc_md_mdkmc_bindunbind.h"

#include "md_bead.h"
#include "md_kmc_bead.h"

void MdMdkmcBindUnbind::Init(space_struct *pSpace, ParticleTracking *pTracking, int ikmc, YAML::Node &node, long seed) {
  KMCBase::Init(pSpace, pTracking, ikmc, node, seed);

  // Grab our specific claims
  eps_eff_  = node["kmc"][ikmc]["eps_eff"].as<double>();
  on_rate_  = node["kmc"][ikmc]["on_rate"].as<double>();
}

void MdMdkmcBindUnbind::Print() {
  printf("MD - MDKMC KMC Module\n");
  KMCBase::Print();
  printf("\t{eps_eff: %2.8f}, {on_rate: %2.8f}\n", eps_eff_, on_rate_);
}

void MdMdkmcBindUnbind::RunKMC(SpeciesBase *spec1, SpeciesBase *spec2) {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();

  // Run the bind unbind
  int g[2] = {0, 1};
  for (int i = 0; i < 2; ++i) {
    int j = gsl_rng_uniform_int(rng_.r, 2);
    int swapme = g[i];
    g[i] = g[j];
    g[j] = swapme;
  }

  if (debug_trace)
    printf("MDKMC module %d -> %d\n", g[0], g[1]);

  for (int i = 0; i < 2; ++i) {
    switch (g[i]) {
      case 0:
        Bind();
        break;
      case 1:
        Unbind(spec1);
        break;
    }
  }

  FinishKMC(spec1);
}

void MdMdkmcBindUnbind::Bind() {
  // Loop over particles and get the appropriate setup
  for (int idx = 0; idx < nsimples_; ++idx) {
    auto part1 = (*simples_)[idx];
    if ((!part1->IsKMC()) || (part1->GetSID() != sid1_)) continue;
    // Dynamic cast to MDKMCBead, we know it has the correct SID
    MDKMCBead *pkmcbead = dynamic_cast<MDKMCBead*>(part1);
    double binding_affinity = eps_eff_ * on_rate_ * pkmcbead->GetDelta();
    auto nexp = pkmcbead->GetNExp();
    if (nexp <  std::numeric_limits<double>::epsilon() &&
        nexp > -std::numeric_limits<double>::epsilon()) nexp = 0.0;
    if (!pkmcbead->GetBound() && nexp > 0.0) {
      auto mrng = pkmcbead->GetRNG();
      double roll = gsl_rng_uniform(mrng->r);
      if (roll < nexp) {
        if (debug_trace)
          printf("[%d] Successful KMC move {bind}, {nexp: %2.4f}, {roll: %2.4f}\n", pkmcbead->GetOID(), nexp, roll);
        pkmcbead->SetBound(true);
        // Figure out where to attach
        double pos = 0.0;
        auto neighbors = tracking_->GetNeighbors();
        for (auto nldx = neighbors[idx].begin(); nldx != neighbors[idx].end(); ++nldx) {
          auto part2 = (*simples_)[nldx->idx_];
          if (part2->GetSID() != sid2_) continue;
          pos += binding_affinity * nldx->kmc_;
          if (pos > roll) {
            // We know this is an MDBead, cast
            MDBead *pbead = dynamic_cast<MDBead*>(part2);
            if (debug_trace)
              printf("[%d,%d] Attaching to [%d,%d] {localpos: %2.4f}\n", idx, pkmcbead->GetOID(), nldx->idx_, pbead->GetOID(), pos);
            pkmcbead->Attach(nldx->idx_);
            break;
          } // found the one to attach to
        } // look @ neighbors
      } // successful kmc move
    } // This one wasn't bound
  } // Loop over all particles, looking for mdkmcbead
}

void MdMdkmcBindUnbind::Unbind(SpeciesBase *spec) {
  double fake_rate = on_rate_ * 10000; // artificially inflate offrate

  // This is done on the species level
  MDKMCBeadSpecies *pkmcbspec = dynamic_cast<MDKMCBeadSpecies*>(spec);
  int nbound = pkmcbspec->GetNBound();
  double poff_single = fake_rate * pkmcbspec->GetDelta();
  int noff = (int)gsl_ran_binomial(rng_.r, poff_single, nbound);
  if (debug_trace)
    printf("[species] {poffsingle: %2.8f, noff: %d}\n", poff_single, noff);
  // Remove noff
  for (int i = 0; i < noff; ++i) {
    int idxloc = -1;
    bool foundidx = false;
    int idxoff = gsl_rng_uniform_int(rng_.r, nbound);
    // Find the one to remove
    for (int idx = 0; idx < nsimples_; ++idx) {
      auto part = (*simples_)[idx];
      if ((!part->IsKMC()) || (part->GetSID() != sid1_)) continue;
      MDKMCBead *pkmcbead = dynamic_cast<MDKMCBead*>(part);
      if (!pkmcbead->GetBound()) continue;
      idxloc++;
      if (idxloc == idxoff) {
        if (debug_trace)
          printf("[%d] Successful KMC move {unbind}, {idxoff=idxloc=%d}\n", pkmcbead->GetOID(), idxloc);
        double randr[3];
        double mag2 = 0.0;
        double mrcut = 0.75;
        double mrcut2 = mrcut * mrcut;
        auto mrng = pkmcbead->GetRNG();
        double prevpos[3];
        std::copy(pkmcbead->GetRigidPosition(), pkmcbead->GetRigidPosition()+ndim_, prevpos);
        do {
          mag2 = 0.0;
          for (int i = 0; i < ndim_; ++i) {
            double mrand = gsl_rng_uniform(mrng->r);
            randr[i] = 2*mrcut*(mrand - 0.5);
            mag2 += SQR(randr[i]);
          }
        } while (mag2 > mrcut2);
        // Randomly set position based on randr
        for (int i = 0; i < ndim_; ++i) {
          randr[i] = randr[i] + prevpos[i];
        }
        pkmcbead->SetPosition(randr);
        pkmcbead->SetBound(false);
        // Set a random velocity
        double newvel[3];
        double newvelpos[3];
        for (int i = 0; i < ndim_; ++i) {
          newvel[i] = 4*(gsl_rng_uniform_pos(mrng->r) -0.5);
          newvelpos[i] = randr[i] - newvel[i]*pkmcbead->GetDelta();
        }
        pkmcbead->SetPrevPosition(newvelpos); 
        if (debug_trace) {
          auto part2 = (*simples_)[pkmcbead->GetAttach()];
          printf("[%d,%d] Detached from [%d,%d] (%2.8f, %2.8f) -> (%2.8f, %2.8f)\n", idx, pkmcbead->GetOID(),
              pkmcbead->GetAttach(), part2->GetOID(),
              part2->GetRigidPosition()[0], part2->GetRigidPosition()[1],
              pkmcbead->GetRigidPosition()[0], pkmcbead->GetRigidPosition()[1]);
        }
        pkmcbead->Attach(-1);
        foundidx = true;
        break;
      } // found the one to detach
    } // loop over all simples

    if (foundidx)
      break;
  } // how many to remove?
}

void MdMdkmcBindUnbind::FinishKMC(SpeciesBase* spec) {
  // Finish the kmc stuff, attach, set positions, etc
  for (int idx = 0; idx < nsimples_; ++idx) {
    auto part = (*simples_)[idx];
    if ((!part->IsKMC()) || (part->GetSID() != sid1_)) continue;
    // Dynamic Cast
    MDKMCBead *pkmcbead = dynamic_cast<MDKMCBead*>(part);
    // If we're bound, update the position to the attached
    if (pkmcbead->GetBound()) {
      pkmcbead->SetNExp(0.0);
      auto aidx = pkmcbead->GetAttach();
      auto part2 = (*simples_)[aidx];
      auto apos = part2->GetRigidPosition();
      auto vpos = part2->GetVelocity();
      auto mpos = pkmcbead->GetRigidPosition();
      if (debug_trace)
        printf("[%d,%d] attached [%d,%d], (%2.4f, %2.4f) -> setting -> (%2.4f, %2.4f)\n",
               idx, pkmcbead->GetOID(), aidx, part2->GetOID(), mpos[0], mpos[1],
               apos[0], apos[1]);
      pkmcbead->SetPrevPosition(mpos);
      pkmcbead->SetPosition(apos);
      pkmcbead->SetVelocity(vpos);
    }
  }

  // XXX CJE add a particle to this species
  //MDKMCBeadSpecies *pkmcbspec = dynamic_cast<MDKMCBeadSpecies*>(spec);
  //pkmcbspec->AddMember();
}
