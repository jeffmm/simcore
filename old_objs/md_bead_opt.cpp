#include "md_bead_opt.h"

//Order must be preserved between Write and read posits
//void MDBeadOpt::WritePosit(std::ofstream &op){
//  //op.write(&oid_, sizeof(oid_), op);
//  for(auto& pos : position_)
//    op.write((char*)&pos, sizeof(pos));
//  for(auto& spos : scaled_position_)
//    op.write((char*)&spos, sizeof(spos));
//  for(auto& u : orientation_)
//    op.write((char*)&u, sizeof(u));
//  op.write((char*)&diameter_, sizeof(diameter_));
//  op.write((char*)&length_, sizeof(length_));
//}

//void MDBeadOpt::ReadPosit(std::ifstream &ip){}

//void MDBeadOptSpecies::WritePosits( std::ofstream &op){
//  for( auto& mem_it : members_){
//    mem_it->WritePosit(op);
//  }
//}

//void MDBeadOptSpecies::ReadPosits(){}

void MDBeadOptSpecies::Configurator() {
  char *filename = params_->config_file;
  std::cout << "MDBeadOpt species\n";

  YAML::Node node = YAML::LoadFile(filename);

  // See what kind of insertion we are doing
  std::string insertion_type;
  insertion_type = node["md_bead_opt"]["properties"]["insertion_type"].as<std::string>();
  std::cout << "   insertion type: " << insertion_type << std::endl;
  bool can_overlap = node["md_bead_opt"]["properties"]["overlap"].as<bool>();
  std::cout << "   overlap: " << (can_overlap ? "true" : "false") << std::endl;

  if (insertion_type.compare("xyz") == 0) {
    std::cout << "Nope, not yet!\n";
    exit(1);
  } else if (insertion_type.compare("random") == 0) {
    int nmdbeads    = node["md_bead_opt"]["mdbead"]["num"].as<int>();
    double diameter = node["md_bead_opt"]["mdbead"]["diameter"].as<double>();
    double mass     = node["md_bead_opt"]["mdbead"]["mass"].as<double>();
    std::cout << "   n_md_beads: " << nmdbeads << std::endl;
    std::cout << "   diameter:   " << diameter << std::endl;
    std::cout << "   mass:       " << mass << std::endl;
    params_->n_md_bead = nmdbeads;
    params_->md_bead_diameter = diameter;
    params_->md_bead_mass = mass;

    for (int i = 0; i < nmdbeads; ++i) {
      MDBeadOpt *member = new MDBeadOpt(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->Init();
      
      // Check if overlaps allowed
      if (can_overlap) {
        members_.push_back(member);
      } else {
        // Check against other md beads
        bool isoverlap = true;
        int numoverlaps = 0;
        do {
          numoverlaps++;
          isoverlap = false;
          for (auto mdit = members_.begin(); mdit != members_.end() && !isoverlap; ++mdit) {
            interactionmindist idm;
            auto part1 = member;
            auto part2 = (*mdit);
            MinimumDistance(part1, part2, idm, space_->n_dim, space_->n_periodic, space_);
            double diameter2 = diameter*diameter;

            if (idm.dr_mag2 < diameter2) {
              isoverlap = true;
              member->Init();
            }
          }
          if (numoverlaps > params_->max_overlap) {
            std::cout << "ERROR: Too many overlaps detected.  Inserted " << i << " of " << nmdbeads;
            std::cout << ".  Check packing ratio for objects.\n";
            exit(1);
          }
        } while (isoverlap);
        members_.push_back(member);
      }
    }

  } else {
    std::cout << "Nope, not yet!\n";
    exit(1);
  }
}
