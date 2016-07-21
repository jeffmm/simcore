// Implementation for simple binding and unbinding

#include "xlink_kmc.h"

#include "xlink.h"
#include "xlink_head.h"
#include "br_rod.h"

// elgacy
#include "br_walker.h"

void XlinkKMC::Init(space_struct *pSpace,
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

void XlinkKMC::Print() {
  printf("Xlink - BR Rod KMC Module\n");
  KMCBase::Print();
  printf("\t{eps_eff: %2.8f}, {on_rate: %2.8f}, {alpha: %2.4f}, {rcut: %2.2f}\n", eps_eff_, on_rate_,
      alpha_, mrcut_);
}

void XlinkKMC::PrepKMC() {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();
  oid_position_map_ = tracking_->GetOIDPositionMap();
  neighbors_ = tracking_->GetNeighbors();

  // Prepare each composite particle for the upcoming kmc step
  if (!spec1_->IsKMC()) return;
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  double ntot_unbound = 0.0;
  auto xlinks = pxspec->GetXlinks(); 

  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    // Determine what state we're in so that we can do the appropriate thing
    // XXX Do the rest of this
    switch((*xit)->GetBoundState()) {
      case unbound:
        Update_0_1(*xit);
        ntot_unbound += (*xit)->GetNExp();
        break;
    }
  }

  pxspec->SetNExp(ntot_unbound);

  return;
}

void XlinkKMC::Update_0_1(Xlink* xit) {
  double nexp_xlink = 0.0;
  auto heads = xit->GetHeads();
  // XXX FIXME
  double binding_affinity = eps_eff_ * on_rate_ * alpha_ * xit->GetDelta();
  for (auto head = heads->begin(); head != heads->end(); ++head) {
    // Look @ neighbors
    double nexp = 0.0;
    auto idx = (*oid_position_map_)[head->GetOID()];
    for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx) {
      nexp += binding_affinity * nldx->kmc_; 
    }
    head->SetNExp(nexp);
    nexp_xlink += nexp;
  }

  xit->SetNExp(nexp_xlink);
}

void XlinkKMC::StepKMC() {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();
  oid_position_map_ = tracking_->GetOIDPositionMap();
  neighbors_ = tracking_->GetNeighbors();

  // Run the bind unbind
  int g[2] = {0, 1};
  for (int i = 0; i < 2; ++i) {
    int j = gsl_rng_uniform_int(rng_.r, 2);
    int swapme = g[i];
    g[i] = g[j];
    g[j] = swapme;
  }

  if (debug_trace)
    printf("XlinkKMC module %d -> %d\n", g[0], g[1]);

  for (int i = 0; i < 2; ++i) {
    switch (g[i]) {
      case 0:
        KMC_0_1();
        break;
      case 1:
        KMC_1_0();
        break;
    }
  }
}

void XlinkKMC::KMC_0_1() {
  // Loop over the xlinks to see who binds
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    // Only take free ones
    if ((*xit)->GetBoundState() != unbound) continue;
    auto nexp = (*xit)->GetNExp();
    if (nexp <  std::numeric_limits<double>::epsilon() &&
        nexp > -std::numeric_limits<double>::epsilon()) nexp = 0.0;
    // IF we have some probability to fall onto a neighbor, check it
    if (nexp > 0.0) {
      auto mrng = (*xit)->GetRNG();
      double roll = gsl_rng_uniform(mrng->r);
      if (roll < nexp) {
        int head_type = gsl_rng_uniform(mrng->r) < 0.5;
        auto heads = (*xit)->GetHeads();
        auto head = heads->begin() + head_type;
        // XXX FIXME
        double binding_affinity = 2.0 * eps_eff_ * on_rate_ * alpha_ * head->GetDelta();
        if (debug_trace)
          printf("[%d] Successful KMC move {0 -> 1}, {nexp: %2.4f}, {roll: %2.4f}, {head: %d}\n", (*xit)->GetOID(),
              nexp, roll, head_type);
        double pos = 0.0;
        int idx = (*oid_position_map_)[head->GetOID()];
        for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx) {
          auto part2 = (*simples_)[nldx->idx_];
          if (part2->GetSID() != sid2_) continue;
          pos += binding_affinity * nldx->kmc_;
          if (pos > roll) {
            if (debug_trace)
              printf("[%d,%d] Attaching to [%d,%d]\n", idx, head->GetOID(), nldx->idx_, part2->GetOID());

            // Here, we do more complicated stuff.  Calculate the coordinate along
            // the rod s.t. the line vector is perpendicular to the separation vec
            // (closest point along carrier line).  In this frame, the position
            // of the xlink should be gaussian distributed
            double r_x[3];
            double r_rod[3];
            double s_rod[3];
            double u_rod[3];
            std::copy(head->GetRigidPosition(), head->GetRigidPosition()+ndim_, r_x);
            std::copy(part2->GetRigidPosition(), part2->GetRigidPosition()+ndim_, r_rod);
            std::copy(part2->GetRigidScaledPosition(), part2->GetRigidScaledPosition()+ndim_, s_rod);
            std::copy(part2->GetRigidOrientation(), part2->GetRigidOrientation()+ndim_, u_rod);
            double l_rod = part2->GetRigidLength();
            double *s_1 = NULL; // FIXME from robert
            double rcontact[3];
            double dr[3];
            double mu = 0.0;
            min_distance_point_carrier_line(ndim_, nperiodic_,
                                            space_->unit_cell, r_x, s_1,
                                            r_rod, s_rod, u_rod, l_rod,
                                            dr, rcontact, &mu);

            double r_min[3];
            double r_min_mag2 = 0.0;
            for (int i = 0; i < ndim_; ++i) {
              r_min[i] = -mu * u_rod[i] - dr[i];
              r_min_mag2 += SQR(r_min[i]);
            }
            double mrcut2 = mrcut_*mrcut_;
            double a = sqrt(mrcut2 - r_min_mag2);
            //double a = sqrt(1.0 - r_min_mag2); //FIXME is this right for 1.0? or mrcut2?
            if (isnan(a))
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
            head->Attach(nldx->idx_, crosspos);
            if (debug_trace)
              printf("\t{mu: %2.4f}, {crosspos: %2.4f}\n", mu, crosspos);
            head->SetBound(true);
            (*xit)->CheckBoundState();
            break;
          } // the actual neighbor to fall on
        } // loop over neighbors to figure out which to attach to
      } // successfully bound
    } // found an xlink with some expectation of binding
  } // loop over all xlinks
}

void XlinkKMC::KMC_1_0() {
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  int nbound1[2];
  std::copy(pxspec->GetNBound1(), pxspec->GetNBound1()+2, nbound1);
  // XXX FIXME fix the on rates
  double poff_single_a = 100 * on_rate_ * alpha_ * pxspec->GetDelta();
  double poff_single_b = 100 * on_rate_ * alpha_ * pxspec->GetDelta();

  int noff[2] = {(int)gsl_ran_binomial(rng_.r, poff_single_a, nbound1[0]),
                 (int)gsl_ran_binomial(rng_.r, poff_single_b, nbound1[1])};
  if (debug_trace)
    printf("[Xlink] {poff_single: (%2.8f, %2.8f)}, {noff: (%d, %d)}\n",
        poff_single_a, poff_single_b, noff[0], noff[1]);

  for (int i = 0; i < (noff[0] + noff[1]); ++i) {
    int head_type = i < noff[1];
    int idxloc = -1;
    int idxoff = gsl_rng_uniform_int(rng_.r, nbound1[head_type]);
    bool foundidx = false;

    // Find the one to remove
    for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
      if ((*xit)->GetBoundState() != singly) continue;
      auto heads = (*xit)->GetHeads();
      // Check to increment only in the case we have the correctly
      // bound head
      if (!(*heads)[head_type].GetBound()) continue;
      idxloc++;
      if (idxloc == idxoff) {
        auto head0 = heads->begin();
        auto head1 = heads->begin()+1;

        // Figure out which head attached
        // Do some fancy aliasing to make this easier
        XlinkHead *attachedhead;
        XlinkHead *nonattachead;
        if (head_type == 0) {
          attachedhead = &(*head0);
          nonattachead = &(*head1);
        } else {
          attachedhead = &(*head1);
          nonattachead = &(*head0);
        }

        if (debug_trace)
          printf("[x:%d,head:%d] Successful KMC move {1 -> 0}, {idxoff=idxloc=%d}, {head: %d}\n",
              (*xit)->GetOID(), attachedhead->GetOID(), idxloc, head_type);
        // Place withint some random distance of the attach point
        double randr[3];
        double mag2 = 0.0;
        double mrcut2 = mrcut_ * mrcut_;
        auto mrng = attachedhead->GetRNG();
        double prevpos[3];
        std::copy(attachedhead->GetRigidPosition(), attachedhead->GetRigidPosition()+ndim_, prevpos);
        do {
          mag2 = 0.0;
          for (int i = 0; i < ndim_; ++i) {
            double mrand = gsl_rng_uniform(mrng->r);
            randr[i] = 2*mrcut_*(mrand - 0.5);
            mag2 += SQR(randr[i]);
          }
        } while(mag2 > mrcut2);
        // Randomly set position based on randr
        for (int i = 0; i < ndim_; ++i) {
          randr[i] = randr[i] + prevpos[i];
        }
        attachedhead->SetPosition(randr);
        attachedhead->SetPrevPosition(prevpos);
        nonattachead->SetPosition(randr);
        nonattachead->SetPrevPosition(prevpos);
        attachedhead->SetBound(false);
        (*xit)->CheckBoundState();

        if (debug_trace) {
          auto attachid = attachedhead->GetAttach();
          auto idx = (*oid_position_map_)[attachedhead->GetOID()];
          auto part2 = (*simples_)[attachid.first];
          printf("[%d,%d] Detached from [%d,%d] (%2.8f, %2.8f) -> (%2.8f, %2.8f)\n",
             idx, attachedhead->GetOID(), attachid.first, part2->GetOID(),
             prevpos[0], prevpos[1],
             attachedhead->GetRigidPosition()[0], attachedhead->GetRigidPosition()[1]);
        }


        attachedhead->Attach(-1, 0.0);
        foundidx = true;
        break;

      } // found it!
    } // find the one to remove

  } // How many to remove
}

void XlinkKMC::UpdateKMC() {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();
  oid_position_map_ = tracking_->GetOIDPositionMap();
  neighbors_ = tracking_->GetNeighbors();

  int nbound_1[2] = {0, 0};
  int nbound_2 = 0;
  int nfree = 0;

  // Do a switch on the type that we're examining
  // Loop over the xlinks to see who  does what
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    switch((*xit)->GetBoundState()) {
      case unbound:
        nfree++;
        break;
      case singly:
        UpdateStage1(*xit, nbound_1);
        break;
    }
  }

  pxspec->SetNFree(nfree);
  pxspec->SetNBound1(nbound_1[0], nbound_1[1]);
  pxspec->SetNBound2(nbound_2);
}

void XlinkKMC::UpdateStage1(Xlink *xit, int *nbound1) {

  // Set nexp to zero for all involved
  xit->SetNExp(0.0);
  auto heads = xit->GetHeads();
  auto head0 = heads->begin();
  auto head1 = heads->begin()+1;
  head0->SetNExp(0.0);
  head1->SetNExp(0.0);
  
  // Figure out which head attached
  // Do some fancy aliasing to make this easier
  XlinkHead *attachedhead;
  XlinkHead *nonattachead;
  if (head0->GetBound()) {
    attachedhead = &(*head0);
    nonattachead = &(*head1);
    nbound1[0]++;
  } else if (head1->GetBound()) {
    attachedhead = &(*head1);
    nonattachead = &(*head0);
    nbound1[1]++;
  } else {
    printf("Something has gone horribly wrong\n");
    exit(1);
  }
  auto aidx = attachedhead->GetAttach().first;
  auto cross_pos = attachedhead->GetAttach().second; // relative to the -end of the rod!!
  auto r_x = attachedhead->GetRigidPosition();
  
  auto part2 = (*simples_)[aidx];
  auto r_rod = part2->GetRigidPosition();
  auto u_rod = part2->GetRigidOrientation();
  auto l_rod = part2->GetRigidLength();

  // If we are moving with some velocity, do that
  cross_pos += velocity_ * attachedhead->GetDelta();
  if (cross_pos > l_rod) {
    cross_pos = l_rod;
  } else if (cross_pos < 0.0) {
    cross_pos = 0.0;
  }
  attachedhead->Attach(aidx, cross_pos);

  double rxnew[3];
  for (int i = 0; i < ndim_; ++i) {
    rxnew[i] = r_rod[i] - 0.5 * u_rod[i] * l_rod + cross_pos * u_rod[i];
  }
  auto idx = (*oid_position_map_)[attachedhead->GetOID()];
  if (debug_trace)
    printf("[%d,%d] attached [%d,%d], (%2.4f, %2.4f) -> setting -> (%2.4f, %2.4f)\n",
           idx, attachedhead->GetOID(), aidx, part2->GetOID(), r_x[0], r_x[1],
           rxnew[0], rxnew[1]);
  attachedhead->SetPrevPosition(r_x);
  attachedhead->SetPosition(rxnew);
  nonattachead->SetPrevPosition(r_x);
  nonattachead->SetPosition(rxnew);
}

void XlinkKMC::Dump() {
  // print out the information appropriate to kmc
  if (debug_trace) {
    XlinkSpecies *pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
    printf("XlinkKMC -> dump\n");
    printf("\t{n_exp(this delta): %2.4f}\n", pxspec->GetNExp());
    printf("\t{nfree:  %d}\n", pxspec->GetNFree());
    printf("\t{nbound1: %d,%d}\n", pxspec->GetNBound1()[0], pxspec->GetNBound1()[1]);
    printf("\t{nbound2: %d}\n", pxspec->GetNBound2());
    pxspec->DumpKMC();
  }
}

void XlinkKMC::PrepOutputs() {
  std::ostringstream file_name;
  file_name << "sc-kmc-XlinkKMC.log";
  kmc_file.open(file_name.str().c_str(), std::ios_base::out);
  kmc_file << "#ntot #nfree #nbound1[0] #nbound1[1] #nbound2\n";
  kmc_file.close();
}

void XlinkKMC::WriteOutputs() {
  XlinkSpecies *pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  std::ostringstream file_name;
  file_name << "sc-kmc-XlinkKMC.log";
  kmc_file.open(file_name.str().c_str(), std::ios_base::out | std::ios_base::app);
  kmc_file.precision(16);
  kmc_file.setf(std::ios::fixed, std::ios::floatfield);
  kmc_file << pxspec->GetNMembers() << " " << pxspec->GetNFree() << " " <<
    pxspec->GetNBound1()[0] << " " << pxspec->GetNBound1()[1] << " " << 
    pxspec->GetNBound2() << "\n";
  kmc_file.close();
}
