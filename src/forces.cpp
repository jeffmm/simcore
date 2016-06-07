#include "forces.h"

void Forces::Init(space_struct *space, std::vector<SpeciesBase*> species, double cell_length) {
  space_=space;
  n_dim_ = space->n_dim;
  n_periodic_ = space->n_periodic;
  cell_list_.Init(n_dim_, n_periodic_, cell_length, space->radius);
  CheckOverlap(species);
  InitPotentials(species);

  // XXX: Cje init my own forces
  test_force_->Init(space, 0.0);
  test_force_->LoadSimples(species);
  test_force_->InitPotentials(species);
}

void Forces::InitPotentials(std::vector<SpeciesBase*> species) {
  for (auto it=species.begin(); it!= species.end(); ++it) {
    auto pot_vec = (*it)->GetPotentials();
    for (auto jt=pot_vec.begin(); jt!=pot_vec.end(); ++jt)
      potentials_.AddPotential(jt->first.first,jt->first.second,jt->second);
  }

  //potentials_.Print();
}

void Forces::DumpAll() {
    // Dump all the particles and their positions, forces, energy (2d)
    #ifdef DEBUG
    for (int i = 0; i < (int)simples_.size(); ++i) {
        auto part = simples_[i];
        auto oid = part->GetOID();
        printf("\to(%d) = ", oid);
        printf("x{%2.2f, %2.2f}, ", part->GetPosition()[0], part->GetPosition()[1]);
        printf("f{%2.2f, %2.2f}, ", part->GetForce()[0], part->GetForce()[1]);
        printf("u{%2.2f}, p{%2.2f}\n", part->GetKineticEnergy(), part->GetPotentialEnergy());
    }
    #endif
}

void Forces::UpdateCellList(std::vector<SpeciesBase*> species) {
  LoadSimples(species);
  interactions_.clear();
  cell_list_.LoadSimples(simples_);
  interactions_ = cell_list_.GetInteractions();
}

void Forces::LoadSimples(std::vector<SpeciesBase*> species) {
  simples_.clear();
  for (auto it=species.begin(); it!=species.end(); ++it) {
    std::vector<Simple*> sim_vec = (*it)->GetSimples();
    simples_.insert(simples_.end(), sim_vec.begin(), sim_vec.end());
  }
}

void Forces::InteractMP() {
    printf("Forces - InteractMP\n");
    test_force_->Interact();
    printf("Forces done - InteractMP\n");
}

void Forces::Interact() {
  for (auto it=interactions_.begin(); it!= interactions_.end(); ++it) {
    auto o1 = it->first;
    auto o2 = it->second;
    if (o1->GetOID() == o2->GetOID()) 
      error_exit("ERROR: Got self-interaction.\n");
    if (o1->GetCID() == o2->GetCID()) continue;
    PotentialBase * pot = potentials_.GetPotential(o1->GetSID(), o2->GetSID());
    if (pot == NULL) continue;
    MinimumDistance(*it);
    if (dr_mag_ > pot->GetRCut()) continue;
    o1->GiveInteraction(FirstInteraction(pot));
    o2->GiveInteraction(SecondInteraction(pot));
  }
}

void Forces::CheckOverlap(std::vector<SpeciesBase*> species) {
  bool overlap = true;
  int num = 0;
  do {
    num++;
    overlap=false;
    UpdateCellList(species);
    for (auto it=interactions_.begin(); it!= interactions_.end(); ++it) {
      if (it->first->GetCID() == it->second->GetCID()) continue;
      MinimumDistance(*it);
      if (dr_mag2_ < buffer_mag2_) {
        overlap = true;
        SID const sid = it->second->GetSID();
        unsigned int const cid = it->second->GetCID();
        for (auto spec=species.begin(); spec!=species.end(); ++spec) {
          if ((*spec)->GetSID() == sid) {
            (*spec)->ReInit(cid);
          }
        }
        break;
      }
    }
    if (num > 10000)
      error_exit("ERROR: Too many overlaps detected. Check packing ratio for objects.\n");
  } while (overlap);
}

void Forces::MinimumDistance(cell_interaction ix) {
  double const * const r1 = ix.first->GetPosition();
  double const * const s1 = ix.first->GetScaledPosition();
  double const * const u1 = ix.first->GetOrientation();
  double const * const r2 = ix.second->GetPosition();
  double const * const s2 = ix.second->GetScaledPosition();
  double const * const u2 = ix.second->GetOrientation();
  double const l1 = ix.first->GetLength();
  double const l2 = ix.second->GetLength();
  double const d1 = ix.first->GetDiameter();
  double const d2 = ix.second->GetDiameter();
  /* TODO: Think about how best to do this for general shapes, like 2d
     polygons that can represent the local surface of more complex 3d
     shapes. Perhaps assume all local surface to be triangular polygons.*/
  dr_mag2_ = 0;
  std::fill(dr_, dr_+3, 0.0);
  std::fill(contact1_, contact1_+3, 0.0);
  std::fill(contact2_, contact2_+3, 0.0);
  buffer_mag_ = 0.5*(d1+d2);
  buffer_mag2_ = buffer_mag_*buffer_mag_;
  if (l1 == 0 && l2 == 0)
    min_distance_point_point(n_dim_, n_periodic_, space_->unit_cell, 
                             r1, s1, r2, s2, dr_, &dr_mag2_);
  else if (l1 == 0 && l2 > 0) 
    min_distance_sphere_sphero(n_dim_, n_periodic_, space_->unit_cell,
                               r1, s1, r2, s2, u2, l2,
                               dr_, &dr_mag2_, contact2_);
  else if (l1 > 0 && l2 == 0) 
    min_distance_sphere_sphero(n_dim_, n_periodic_, space_->unit_cell,
                               r2, s2, r1, s1, u1, l1,
                               dr_, &dr_mag2_, contact1_);
  else if (l1 > 0 && l2 > 0)
    min_distance_sphero(n_dim_, n_periodic_, space_->unit_cell,
                        r1, s1, u1, l1, r2, s2, u2, l2,
                        dr_, &dr_mag2_, contact1_, contact2_);
  dr_mag_ = sqrt(dr_mag2_);
}

interaction Forces::FirstInteraction(PotentialBase *pot) {
  interaction ix;
  std::copy(dr_, dr_+3, ix.dr);
  std::copy(contact1_, contact1_+3, ix.contact);
  ix.buffer = buffer_mag_;
  ix.dr_mag = dr_mag_;
  ix.potential = pot;
  return ix;
}

interaction Forces::SecondInteraction(PotentialBase *pot) {
  interaction ix;
  std::copy(dr_, dr_+3, ix.dr);
  // switch direction of dr
  for (int i=0; i<3; ++i)
    ix.dr[i] = -ix.dr[i];
  std::copy(contact2_, contact2_+3, ix.contact);
  ix.buffer = buffer_mag_;
  ix.dr_mag = dr_mag_;
  ix.potential = pot;
  return ix;
}

