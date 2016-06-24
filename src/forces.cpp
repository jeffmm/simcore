#include <chrono>

#include "forces.h"
#include "object.h"

// Pass in the main system properties information
void Forces::Init(system_parameters *pParams, space_struct *pSpace, std::vector<SpeciesBase*> species) {
  params_ = pParams;
  space_ = pSpace;
  n_dim_ = space_->n_dim;
  n_periodic_ = space_->n_periodic;
  draw_flag_ = params_->draw_interactions;
  //XXX: CJE switch on force type for now
  printf("********\n");
  printf("Create Forces\n");
  switch (params_->ftype) {
    case 0:
        printf("Must specify a force substructure, exiting!\n");
        exit(1);
    case 1:
        printf("->Using all pairs (brute) force substructure\n");
        force_type_ = FTYPE::allpairs;
        force_module_ = forceFactory<ForceBrute>();
        skin_ = 0.0;
        break;
    case 2:
        printf("->Using microcells force substructure\n");
        force_type_ = FTYPE::microcells;
        force_module_ = forceFactory<ForceMicrocell>();
        skin_ = 0.0;
        printf("ERROR: Currently deprecated on account of composite objects, exiting\n");
        exit(1);
        break;
    case 3:
        printf("->Using cells force substructure\n");
        force_type_ = FTYPE::cells;
        force_module_ = forceFactory<ForceCell>();
        skin_ = 0.0;
        printf("ERROR: Currently deprecated on account of composite objects, exiting\n");
        exit(1);
        break;
    case 4:
        printf("->Using neighbor lists all pairs substructure\n");
        force_type_ = FTYPE::neighborallpairs;
        force_module_ = forceFactory<ForceNeighborListAP>();
        skin_ = params_->masterskin;
        break;
    case 5:
        printf("->Using neighbor list cells substructure\n");
        force_type_ = FTYPE::neighborcells;
        force_module_ = forceFactory<ForceNeighborListCells>();
        skin_ = params_->masterskin;
        printf("ERROR: Currently deprecated on account of composite objects, exiting\n");
        exit(1);
        break;
    default:
        printf("Must specify a force substructure, exiting!\n");
        break;
   }
  //XXX: CJE needed for the check overlaps, change this
  //cell_list_.Init(n_dim_, n_periodic_, params_->cell_length, space_->radius);
  InitPotentials(species);
  CheckOverlap(species);

  // Create the force submodule to do the calculations!
  force_module_->Init(space_, skin_);
  force_module_->LoadSimples(species);
  force_module_->InitPotentials(&potentials_);
  force_module_->Finalize();
  force_module_->UpdateScheme();
  force_module_->print();
}

void Forces::InitPotentials(std::vector<SpeciesBase*> species) {
  // Ask the potential manager to parse the potentials file
  potentials_.Init(space_, params_->potfile);
}

void Forces::DumpAll() {
    // Dump all the particles and their positions, forces, energy (2d)
    #ifdef DEBUG
    force_module_->dump();
    #endif
}

void Forces::UpdateCellList(std::vector<SpeciesBase*> species) {
  LoadSimples(species);
  interactions_.clear();
  cell_list_.LoadSimples(simples_);
  interactions_ = cell_list_.GetInteractions();
}

void Forces::UpdateScheme(std::vector<SpeciesBase*> species) {
  LoadSimples(species);
  force_module_->UpdateScheme();
}

void Forces::LoadSimples(std::vector<SpeciesBase*> species) {
  simples_.clear();
  for (auto it=species.begin(); it!=species.end(); ++it) {
    std::vector<Simple*> sim_vec = (*it)->GetSimples();
    simples_.insert(simples_.end(), sim_vec.begin(), sim_vec.end());
  }
}

void Forces::InteractMP() {
    force_module_->Interact();
}

void Forces::Interact() {
  if (draw_flag_) {
    draw_ = false;
    draw_array_.clear();
  }
  for (auto it=interactions_.begin(); it!= interactions_.end(); ++it) {
    auto o1 = it->first;
    auto o2 = it->second;
    if (o1->GetOID() == o2->GetOID()) 
      error_exit("ERROR: Got self-interaction.\n");
    if (o1->GetCID() == o2->GetCID()) continue;
    PotentialBase * pot = potentials_.GetPotential(o1->GetSID(), o2->GetSID());
    if (pot == NULL) continue;
    //MinimumDistance(*it);
    if (dr_mag_ > pot->GetRCut()) continue;
    o1->GiveInteraction(FirstInteraction(pot));
    o2->GiveInteraction(SecondInteraction(pot));
  }
}

void Forces::CheckOverlap(std::vector<SpeciesBase*> species) {
    // very yucky brute force method of checking overlaps
    // load all the particles

    bool overlap = true;
    int num = 0;
    do {
        num++;
        overlap = false;
        LoadSimples(species);
        for (int idx = 0; idx < simples_.size() - 1 && !overlap; ++idx) {
            for (int jdx = idx + 1; jdx < simples_.size() && !overlap; ++jdx) {
                auto part1 = simples_[idx];
                auto part2 = simples_[jdx];

                if (part1->GetCID() == part2->GetCID()) continue;

                // Check to see if interacting
                PotentialBase *pot = potentials_.GetPotential(part1->GetSID(), part2->GetSID());
                if (pot == nullptr) continue;

                // Check to see if overlap on this is allowed (KMC stuff)
                if (pot->CanOverlap()) continue;

                // Can possibly overlap, check
                interactionmindist idm;
                MinimumDistance(part1, part2, idm, n_dim_, n_periodic_, space_);
                if (idm.dr_mag2 < idm.buffer_mag2) {
                    if (debug_trace) {
                        printf("Overlap detected [oid: %d,%d], [idx: %d, %d], [sid: %d, %d], [cid: %d, %d] -> (%2.2f, %2.2f)\n",
                                part1->GetOID(), part2->GetOID(), idx, jdx, part1->GetSID(), part2->GetSID(), 
                                part1->GetCID(), part2->GetCID(), idm.buffer_mag2, idm.dr_mag2);
                    }
                    overlap = true;
                    auto sid = part2->GetSID();
                    auto cid = part2->GetCID();
                    for (auto spec = species.begin(); spec != species.end(); ++spec) {
                        if ((*spec)->GetSID() == sid) {
                            (*spec)->ReInit(cid);
                        }
                    }
                }
            }
        }
        if (num > 10000)
            error_exit("ERROR: Too many overlaps detected. Check packing ratio for objects.\n");
    } while (overlap);
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

void Forces::Draw(std::vector<graph_struct*> * graph_array) {
  if (!draw_)
    return;
  for (auto it=draw_array_.begin(); it!=draw_array_.end(); ++it)
    graph_array->push_back(&(*it));
}


