// Implementation for simple binding and unbinding

#include "kmc_br_bindunbind.h"

#include "br_walker.h"
#include "br_rod.h"

void BrBindUnbind::Init(space_struct *pSpace,
                        ParticleTracking *pTracking,
                        SpeciesBase *spec1,
                        SpeciesBase *spec2,
                        int ikmc,
                        YAML::Node &node,
                        long seed) {
  KMCBase::Init(pSpace, pTracking, spec1, spec2, ikmc, node, seed);

  // Grab our specific claims
  eps_eff_  = node["kmc"][ikmc]["eps_eff"].as<double>();
  on_rate_  = node["kmc"][ikmc]["on_rate"].as<double>();
  alpha_    = node["kmc"][ikmc]["alpha"].as<double>();
  mrcut_    = node["kmc"][ikmc]["rcut"].as<double>();
  velocity_ = node["kmc"][ikmc]["velocity"].as<double>();
}

void BrBindUnbind::Print() {
  printf("BR Walker - BR Rod KMC Module\n");
  KMCBase::Print();
  printf("\t{eps_eff: %2.8f}, {on_rate: %2.8f}, {alpha: %2.4f}, {rcut: %2.2f}\n", eps_eff_, on_rate_,
      alpha_, mrcut_);
}

void BrBindUnbind::PrepKMC() {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();

  // Prepare each particle/species for the upcoming kmc step
  if (!spec1_->IsKMC()) return;
  BrWalkerSpecies* pwspec = dynamic_cast<BrWalkerSpecies*>(spec1_);
  double ntot = 0.0;
  for (int idx = 0; idx < nsimples_; ++idx) {
    auto part1 = (*simples_)[idx];
    if (!part1->IsKMC() || (part1->GetSID() != sid1_)) continue;
    // We know we have a BrWalker
    BrWalker *pwalker = dynamic_cast<BrWalker*>(part1);
    pwalker->SetNExp(0.0);
    if (pwalker->GetBound()) continue;
    auto neighbors = tracking_->GetNeighbors();
    double binding_affinity = eps_eff_ * on_rate_ * alpha_ * pwalker->GetDelta();
    double nexp = 0.0;
    for (auto nldx = neighbors[idx].begin(); nldx != neighbors[idx].end(); ++nldx) {
      nexp += binding_affinity * nldx->kmc_;
    }
    pwalker->SetNExp(nexp);
    ntot += nexp;
  }

  pwspec->SetNExp(ntot);
}

void BrBindUnbind::StepKMC() {
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
    printf("BR ROD module %d -> %d\n", g[0], g[1]);

  for (int i = 0; i < 2; ++i) {
    switch (g[i]) {
      case 0:
        Bind();
        break;
      case 1:
        Unbind();
        break;
    }
  }
}

void BrBindUnbind::Bind() {
  // Loop over particles and get the appropriate setup
  for (int idx = 0; idx < nsimples_; ++idx) {
    auto part1 = (*simples_)[idx];
    if ((!part1->IsKMC()) || (part1->GetSID() != sid1_)) continue;
    // Dynamic cast to MDKMCBead, we know it has the correct SID
    BrWalker *pwalker = dynamic_cast<BrWalker*>(part1);
    double binding_affinity = eps_eff_ * on_rate_ * alpha_ * pwalker->GetDelta();
    auto nexp = pwalker->GetNExp();
    if (nexp <  std::numeric_limits<double>::epsilon() &&
        nexp > -std::numeric_limits<double>::epsilon()) nexp = 0.0;
    if (!pwalker->GetBound() && nexp > 0.0) {
      auto mrng = pwalker->GetRNG();
      double roll = gsl_rng_uniform(mrng->r);
      if (roll < nexp) {
        if (debug_trace)
          printf("[%d] Successful KMC move {bind}, {nexp: %2.4f}, {roll: %2.4f}\n", pwalker->GetOID(), nexp, roll);
        pwalker->SetBound(true);
        // Figure out where to attach
        double pos = 0.0;
        auto neighbors = tracking_->GetNeighbors();
        for (auto nldx = neighbors[idx].begin(); nldx != neighbors[idx].end(); ++nldx) {
          auto part2 = (*simples_)[nldx->idx_];
          if (part2->GetSID() != sid2_) continue;
          pos += binding_affinity * nldx->kmc_;
          if (pos > roll) {
            // This must remain a simple (groan)
            if (debug_trace)
              printf("[%d,%d] Attaching to [%d,%d] {localpos: %2.4f}\n", idx, pwalker->GetOID(), nldx->idx_, part2->GetOID(), pos);

            // Here, we do more complicated stuff.  Calculate the coordinate along
            // the rod s.t. the line vector is perpendicular to the separation vec
            // (closest point along carrier line).  In this frame, the position
            // of the xlink should be gaussian distributed
            double r_x[3];
            double r_rod[3];
            double s_rod[3];
            double u_rod[3];
            std::copy(pwalker->GetRigidPosition(), pwalker->GetRigidPosition()+ndim_, r_x);
            std::copy(part2->GetRigidPosition(), part2->GetRigidPosition()+ndim_, r_rod);
            std::copy(part2->GetRigidScaledPosition(), part2->GetRigidScaledPosition()+ndim_, s_rod);
            std::copy(part2->GetRigidOrientation(), part2->GetRigidOrientation()+ndim_, u_rod);
            double l_rod = part2->GetRigidLength();
            double *s_1 = NULL; // FIXME from robert
            double rcontact[3];
            double dr[3];
            min_distance_point_carrier_line(ndim_, nperiodic_,
                                            space_->unit_cell, r_x, s_1,
                                            r_rod, s_rod, u_rod, l_rod,
                                            dr, rcontact);

            // Back calculate mu
            double mu = 0.0;
            for (int i = 0; i < ndim_; ++i) {
              if (u_rod[i] > 0.0) {
                mu = rcontact[i]/u_rod[i];
                break;
              }
            }
            double r_min[3];
            double r_min_mag2 = 0.0;
            for (int i = 0; i < ndim_; ++i) {
              r_min[i] = -mu * u_rod[i] - dr[i];
              r_min_mag2 += SQR(r_min[i]);
            }
            double mrcut2 = mrcut_*mrcut_;
            double a = sqrt(1.0 - r_min_mag2); //FIXME is this right for 1.0? or mrcut2?
            if (std::isnan(a))
              a = 0.0;

            double crosspos = 0.0;
            for (int i = 0; i < 100; ++i) {
              double uroll = gsl_rng_uniform(mrng->r);
              crosspos = (uroll - 0.5)*a + mu + 0.5*l_rod;

              if (crosspos >= 0 && crosspos <= l_rod)
                break;
              if (i == 99) {
                crosspos = -mu + 0.5*l_rod;
                if (crosspos < 0)
                  crosspos = 0;
                else if (crosspos > l_rod)
                  crosspos = l_rod;
              }
            }
            pwalker->Attach(nldx->idx_, crosspos); 
            break;
          } // found the one to attach to
        } // look @ neighbors
      } // successful kmc move
    } // This one wasn't bound
  } // Loop over all particles, looking for mdkmcbead
}

void BrBindUnbind::Unbind() {
  // This is done on the species level
  BrWalkerSpecies *pwspec = dynamic_cast<BrWalkerSpecies*>(spec1_);
  int nbound = pwspec->GetNBound();
  double poff_single = on_rate_ * alpha_ * pwspec->GetDelta();
  int noff = (int)gsl_ran_binomial(rng_.r, poff_single, nbound);
  if (debug_trace)
    printf("[BrWalkerSpecies] {poffsingle: %2.8f, noff: %d}\n", poff_single, noff);
  // Remove noff
  for (int i = 0; i < noff; ++i) {
    int idxloc = -1;
    bool foundidx = false;
    int idxoff = gsl_rng_uniform_int(rng_.r, nbound);
    // Find the one to remove
    for (int idx = 0; idx < nsimples_; ++idx) {
      auto part = (*simples_)[idx];
      if ((!part->IsKMC()) || (part->GetSID() != sid1_)) continue;
      BrWalker *pwalker = dynamic_cast<BrWalker*>(part);
      if (!pwalker->GetBound()) continue;
      idxloc++;
      if (idxloc == idxoff) {
        if (debug_trace)
          printf("[%d] Successful KMC move {unbind}, {idxoff=idxloc=%d}\n", pwalker->GetOID(), idxloc);
        double randr[3];
        double mag2 = 0.0;
        double mrcut2 = mrcut_ * mrcut_;
        auto mrng = pwalker->GetRNG();
        double prevpos[3];
        std::copy(pwalker->GetRigidPosition(), pwalker->GetRigidPosition()+ndim_, prevpos);
        do {
          mag2 = 0.0;
          for (int i = 0; i < ndim_; ++i) {
            double mrand = gsl_rng_uniform(mrng->r);
            randr[i] = 2*mrcut_*(mrand - 0.5);
            mag2 += SQR(randr[i]);
          }
        } while (mag2 > mrcut2);
        // Randomly set position based on randr
        for (int i = 0; i < ndim_; ++i) {
          randr[i] = randr[i] + prevpos[i];
        }
        pwalker->SetPosition(randr);
        pwalker->SetPrevPosition(prevpos);
        pwalker->SetBound(false);
        if (debug_trace) {
          auto attachid = pwalker->GetAttach();
          auto part2 = (*simples_)[attachid.first];
          printf("[%d,%d] Detached from [%d,%d] (%2.8f, %2.8f) -> (%2.8f, %2.8f)\n", idx, pwalker->GetOID(),
              attachid.first, part2->GetOID(),
              part2->GetRigidPosition()[0], part2->GetRigidPosition()[1],
              pwalker->GetRigidPosition()[0], pwalker->GetRigidPosition()[1]);
        }
        pwalker->Attach(-1, 0.0);
        foundidx = true;
        break;
      } // found the one to detach
    } // loop over all simples to detach

  } // how many to remove?
}

void BrBindUnbind::UpdateKMC() {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();

  int nbound = 0;
  int nfree = 0;

  // Finish the kmc stuff, attach, set positions, etc
  for (int idx = 0; idx < nsimples_; ++idx) {
    auto part = (*simples_)[idx];
    if ((!part->IsKMC()) || (part->GetSID() != sid1_)) continue;
    // Dynamic Cast
    BrWalker *pwalker = dynamic_cast<BrWalker*>(part);
    // If we're bound, update the position to the attached
    if (pwalker->GetBound()) {
      pwalker->SetNExp(0.0);
      auto aidx = pwalker->GetAttach().first;
      auto cross_pos = pwalker->GetAttach().second; // relative to the -end of the rod!!
      auto r_x = pwalker->GetRigidPosition();
      
      auto part2 = (*simples_)[aidx];
      auto r_rod = part2->GetRigidPosition();
      auto u_rod = part2->GetRigidOrientation();
      auto l_rod = part2->GetRigidLength();

      // If we are moving with some velocity, do that
      cross_pos += velocity_ * pwalker->GetDelta();
      if (cross_pos > l_rod) {
        cross_pos = l_rod;
      } else if (cross_pos < 0.0) {
        cross_pos = 0.0;
      }
      pwalker->Attach(aidx, cross_pos);

      double rxnew[3];
      for (int i = 0; i < ndim_; ++i) {
        rxnew[i] = r_rod[i] - 0.5 * u_rod[i] * l_rod + cross_pos * u_rod[i];
      }

      if (debug_trace)
        printf("[%d,%d] attached [%d,%d], (%2.4f, %2.4f) -> setting -> (%2.4f, %2.4f)\n",
               idx, pwalker->GetOID(), aidx, part2->GetOID(), r_x[0], r_x[1],
               rxnew[0], rxnew[1]);
      pwalker->SetPrevPosition(r_x);
      pwalker->SetPosition(rxnew);
      nbound++;
    } else {
      nfree++;
    }
  }

  BrWalkerSpecies *pwspec = dynamic_cast<BrWalkerSpecies*>(spec1_);
  pwspec->SetNBound(nbound);
  pwspec->SetNFree(nfree);
}

void BrBindUnbind::Dump() {
  // print out the information appropriate to kmc
  if (debug_trace) {
    BrWalkerSpecies *pwspec = dynamic_cast<BrWalkerSpecies*>(spec1_);
    printf("BrBindUnbind -> dump\n");
    printf("\t{n_exp(this delta): %2.4f}\n", pwspec->GetNExp());
    printf("\t{nfree:  %d}\n", pwspec->GetNFree());
    printf("\t{nbound: %d}\n", pwspec->GetNBound());
    pwspec->DumpKMC();
  }
}
