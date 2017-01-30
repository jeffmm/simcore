#include "br_bead.h"

void BrBead::Init() {
  Simple::Init();
  std::fill(body_frame_,body_frame_+6,0);
}

void BrBead::UpdatePosition() {
  Translate();
  Rotate();
  UpdatePeriodic();
}

void BrBead::KickBead() {
  for (int i=0; i<n_dim_; ++i) {
    double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
    force_[i] += kick*diffusion_;
  }
}

void BrBead::Translate() {
  // Add random component to force
  KickBead();
  Driving();
  std::copy(position_, position_+n_dim_, prev_position_);
  for (int i = 0; i < n_dim_; ++i) {
    position_[i] = position_[i] + force_[i] * delta_ / diameter_;
    dr_tot_[i] += position_[i] - prev_position_[i];
  }
}

void BrBead::Driving() {
  double f_dr[3];
  for (int i=0; i<n_dim_; ++i)
    f_dr[i] = orientation_[i]*driving_factor_;
  AddForce(f_dr);
}

void BrBead::SetDiffusion() {
  friction_rot_ = 3.0*diameter_*diameter_*diameter_;
  rand_sigma_rot_ = sqrt(2.0*delta_/friction_rot_);
  diffusion_ = sqrt(24.0*diameter_/delta_);
}

void BrBead::Rotate() {
// Now handle the random orientation update
  GetBodyFrame();
  for (int j=0; j<n_dim_-1; ++j) {
    double mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_rot_);
    for (int i=0; i<n_dim_; ++i)
      orientation_[i] += mag * body_frame_[n_dim_*j+i];
  }
  normalize_vector(orientation_, n_dim_);
}

/* calculates vector(s) orthogonal to orientation */
void BrBead::GetBodyFrame() {
  if (n_dim_==2) {
    body_frame_[0] = orientation_[1];
    body_frame_[1] = -orientation_[0];
  }
  else {
    double vect1[3] = {1.0, 0.0, 0.0};
    double vect2[3] = {0.0, 1.0, 0.0};
    if (1.0 - ABS(orientation_[0]) > 1e-2)
      cross_product(orientation_, vect1, &(body_frame_[0]), n_dim_);
    else
      cross_product(orientation_, vect2, &(body_frame_[0]), n_dim_);
    normalize_vector(&(body_frame_[0]),n_dim_);
    cross_product(orientation_, &(body_frame_[0]), &(body_frame_[3]), n_dim_);
  }
}

void BrBeadSpecies::Configurator() {
  std::cout << "BRBead species\n";

  // See what kind of insertion we are doing
  std::string insertion_type;
  double color[4] = {1.0, 0.0, 0.0, 1.0};
  insertion_type = params_->br_bead.insertion_type;
  std::cout << "   insertion type: " << insertion_type << std::endl;
  bool can_overlap = params_->br_bead.overlap;
  std::cout << "   overlap: " << (can_overlap ? "true" : "false") << std::endl;
  for (int i = 0; i < 4; ++i) {
    color[i] = params_->br_bead.color;
  }
  std::cout << "   color: [" << color[0] << ", " << color[1] << ", " << color[2] << ", "
    << color[3] << "]\n";
  int draw_type = 0;
  std::string draw_type_s = params_->br_bead.draw_type;
  std::cout << "   draw_type: " << draw_type_s << std::endl;
  if (draw_type_s.compare("flat") == 0) {
    draw_type = 0;
  } else if (draw_type_s.compare("orientation") == 0) {
    draw_type = 1;
  } else {
    draw_type = 2;
  }

  if (insertion_type.compare("xyz") == 0) {
    std::cout << "Nope, not yet!\n";
    exit(1);
  } else if (insertion_type.compare("random") == 0) {

    int n_beads = params_->br_bead.num;
    double diameter = params_->br_bead.diameter;
    std::cout << "   n_br_beads: " << n_beads << std::endl;
    std::cout << "   diameter:   " << diameter << std::endl;

    for (int i = 0; i < n_beads; ++i) {
      BrBead *member = new BrBead(params_, space_, gsl_rng_get(rng_.r), GetSID());
      member->Init();
      //member->SetColor(color, draw_type);
      
      // Check if overlaps allowed
      if (can_overlap) {
        members_.push_back(member);
      } else {
        // Check against other br beads
        bool isoverlap = true;
        int numoverlaps = 0;
        do {
          numoverlaps++;
          isoverlap = false;
          for (auto brit = members_.begin(); brit != members_.end() && !isoverlap; ++brit) {
            Interaction idm;
            auto part1 = member;
            auto part2 = (*brit);
            MinimumDistance(part1, part2, &idm, space_);
            double diameter2 = diameter*diameter;

            if (idm.dr_mag2 < diameter2) {
              isoverlap = true;
              member->Init();
            }
          }
          if (numoverlaps > params_->max_overlap) {
            std::cout << "ERROR: Too many overlaps detected.  Inserted " << i << " of " << n_beads;
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

