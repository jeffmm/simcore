#include "spindle.h"

Spindle::Spindle() : Object() {
  color_ = params_->spindle.color;
  draw_ = draw_type::_from_string(params_->spindle.draw_type.c_str());
  diameter_ = params_->spindle.diameter;
  n_filaments_ = n_filaments_min_ = n_filaments_max_ = params_->spindle.n_filaments;
  k_spring_ = params_->spindle.k_spring;
  k_align_ = params_->spindle.k_align;
  spring_length_ = params_->spindle.spring_length;
  alignment_potential_ = (params_->spindle.alignment_potential ? true : false);
  fixed_spacing_ = false;
  anchor_distance_ = 0.5*diameter_;
  spb_radius_ = params_->spindle.spb_radius;
  SetDiffusion();
}

void Spindle::Init() {
  if (n_filaments_min_ > n_filaments_max_) {
    warning("Min n_filaments greater than max n_filaments in spindle. Setting equal.");
  }
  filaments_.reserve(n_filaments_);
  anchors_.resize(n_filaments_);
  bool out_of_bounds;
  int n_insert = 0;
  do {
    InsertSpindle();
    double buffer = MIN(1.5*params_->filament.min_length,1.1*params_->filament.max_length);
    if (out_of_bounds = CheckBounds(buffer)) continue;
    GenerateAnchorSites();
    filaments_.clear();
    for (int i=0;i<n_filaments_;++i) {
      n_insert = 0;
      do {
        Filament fil;
        filaments_.push_back(fil);
        filaments_.back().SetAnchor(&anchors_[i]);
        filaments_.back().Init(true);
        if (out_of_bounds = filaments_.back().CheckBounds()) {
          filaments_.pop_back();
        }
        if (out_of_bounds && n_insert++ > 100) {
          //if (!fixed_spacing_) RandomizeAnchorPosition(i);
          if (n_insert > 1000) {
            break;
          }
        }
      } while (out_of_bounds);
      if (out_of_bounds) {
        break;
      }
    }
  } while (out_of_bounds);
}

int Spindle::GetCount() {
  int count=1;
  for (filament_iterator it=filaments_.begin();it!=filaments_.end();++it) {
    count += it->GetCount();
  }
  return count;
}

void Spindle::ZeroForce() {
  Object::ZeroForce();
  for (filament_iterator it=filaments_.begin();it!=filaments_.end();++it) {
    it->ZeroForce();
  }
  for (anchor_iterator it=anchors_.begin(); it!=anchors_.end(); ++it) {
    it->ZeroForce();
  }
}

void Spindle::InsertSpindle() {
  if (params_->spindle.insertion_type.compare("random") == 0) {
    InsertRandom();
  }
  else if (params_->spindle.insertion_type.compare("random_oriented") == 0) {
    InsertRandom();
    std::fill(orientation_,orientation_+3,0.0);
    orientation_[n_dim_-1] = 1.0;
  }
  else if (params_->spindle.insertion_type.compare("centered_random") == 0) {
    std::fill(position_,position_+3,0.0);
    generate_random_unit_vector(n_dim_,orientation_,rng_.r);
  }
  else if (params_->spindle.insertion_type.compare("centered_oriented") == 0 ) {
    std::fill(position_,position_+3,0.0);
    std::fill(orientation_,orientation_+3,0.0);
    orientation_[n_dim_-1] = 1.0;
  }
  else {
    error_exit("Spindle insertion type not recognized!");
  }
}

void Spindle::GenerateAnchorSites() {
  anchor_distance_ += (0.5 + 1e-3)*params_->filament.diameter;
  double theta = 0.0;
  double dtheta = 2.0*M_PI/(n_filaments_);
  for (int i_fil=0;i_fil<n_filaments_;++i_fil) {
    if (fixed_spacing_ && n_dim_ == 2) {
      anchors_[i_fil].orientation_[0] = cos(theta);
      anchors_[i_fil].orientation_[1] = sin(theta);
      theta+=dtheta;
    }
    else if (fixed_spacing_ && n_dim_ == 3) {
      warning("Fixed filament spacing not yet implemented for 3D in spindle. Inserting randomly.");
      generate_random_unit_vector(n_dim_,anchors_[i_fil].orientation_,rng_.r);
    }
    else {
      generate_random_unit_vector(n_dim_,anchors_[i_fil].orientation_,rng_.r);
    }
    for (int i=0;i<n_dim_;++i) {
      anchors_[i_fil].position_[i] = position_[i] + anchor_distance_*anchors_[i_fil].orientation_[i];
    }
    anchors_[i_fil].k_spring_ = k_spring_;
    anchors_[i_fil].k_align_ = k_align_;
    anchors_[i_fil].spring_length_ = spring_length_;
    anchors_[i_fil].alignment_potential_ = alignment_potential_;
  }
}

void Spindle::RandomizeAnchorPosition(int i_fil) {
  generate_random_unit_vector(n_dim_,anchors_[i_fil].orientation_,rng_.r);
  for (int i=0;i<n_dim_;++i) {
    anchors_[i_fil].position_[i] = position_[i] + anchor_distance_*anchors_[i_fil].orientation_[i];
  }
}

void Spindle::UpdatePosition(bool midstep) {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<std::vector<Filament>::iterator, std::vector<Filament>::iterator> > chunks;
  chunks.reserve(max_threads); 
  size_t chunk_size= filaments_.size() / max_threads;
  filament_iterator cur_iter = filaments_.begin();
  for(int i = 0; i < max_threads - 1; ++i) {
    filament_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, filaments_.end()));
#pragma omp parallel shared(chunks)
  {
#pragma omp for 
    for(int i = 0; i < max_threads; ++i) {
      for(auto it = chunks[i].first; it != chunks[i].second; ++it) {
        it->UpdatePosition(midstep);
      }
    }
  }
#else
  for (filament_iterator it=filaments_.begin(); it!=filaments_.end(); ++it) {
    it->UpdatePosition(midstep);
  }
#endif
  if (!midstep) {
    SetPrevPosition(position_);
    ApplyForcesTorques();
    Integrate();
    UpdatePeriodic();
  }
}

void Spindle::ApplyForcesTorques() {
  for (int i=0; i<n_dim_; ++i) {
    double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
    force_[i] += kick*diffusion_;
  }
  for (anchor_iterator it=anchors_.begin(); it!=anchors_.end(); ++it) {
    // First calculate translational forces
    //double f_mag = dot_product(n_dim_, it->orientation_, it->force_);
    double dr[3] = {0,0,0};
    for (int i=0; i<n_dim_; ++i) {
      //force_[i] += f_mag * it->orientation_[i];
      force_[i] += it->force_[i];
      dr[i] = it->orientation_[i]*anchor_distance_;
    }
    // Then calculate torques
    double i_torque[3] = {0,0,0};
    cross_product(dr, it->force_, i_torque, n_dim_);
    // 3-dimensions due to torques
    for (int i=0; i<3; ++i) {
      torque_[i] += i_torque[i];
      // And add torque from alignment potential
      if (alignment_potential_) {
        torque_[i] += it->torque_[i];
      }
    }
  }
}

void Spindle::SetDiffusion() {
  gamma_trans_ = 1.0/(diameter_);
  gamma_rot_ = 3.0/CUBE(diameter_);
  diffusion_ = sqrt(24.0*diameter_/delta_);
}

void Spindle::Translate() {
  //double f_cutoff = 0.1 / gamma_trans_ / delta; // ????? arbitrary holdovers from snowman
  //double force_mag = 0.0;
  //for (int i=0; i<n_dim_; ++i)
    //force_mag += SQR(force_[i]);
  //if (force_mag > SQR(f_cutoff)) {
    //normalize_vector(force_,n_dim);
    //for (int i=0; i<n_dim_; ++i)
      //force_[i]*=f_cutoff;
  //}
  double dr[3];
  for (int i=0; i<n_dim_; ++i) {
    dr[i] = force_[i] * delta_ * gamma_trans_;
    position_[i] += dr[i];
  }
  for (anchor_iterator it=anchors_.begin(); it!=anchors_.end(); ++it) {
    for (int i=0; i<n_dim_; ++i)
      it->position_[i]+=dr[i];
  }
}
void Spindle::Rotate() {
  double unit_torque[3], temp[3], r_rel[3];
  double domega, cos_domega, sin_domega, torque_mag;
  // First rotate orientation vector of sphere
  if (n_dim_ == 2) {
    domega = torque_[2]*delta_*gamma_rot_;
    cos_domega = cos(domega);
    sin_domega = sin(domega);
    std::copy(orientation_,orientation_+3,temp);
    orientation_[0] = cos_domega * temp[0] - sin_domega * temp[1];
    orientation_[1] = sin_domega * temp[0] + cos_domega * temp[1];
  }
  else if (n_dim_ == 3) {
    torque_mag = 0.0;
    for (int i=0; i<3; ++i) 
      torque_mag += torque_[i];
    for (int i=0; i<3; ++i) 
      unit_torque[i] = torque_[i]/torque_mag;
    domega = torque_mag * delta_*gamma_rot_;
    rotate_3d_vector(domega, orientation_, unit_torque);
  }
  normalize_vector(orientation_, n_dim_);
  // Now rotate all the attachment sites of the attached filaments
  for (anchor_iterator it=anchors_.begin(); it!=anchors_.end(); ++it) {
    for (int i=0; i<n_dim_; ++i)
      r_rel[i]  = it->position_[i] - position_[i];
    if (n_dim_ == 2) {
      for (int i=0; i<n_dim_; ++i) 
        temp[i] = r_rel[i];
      r_rel[0] = cos_domega * temp[0] - sin_domega * temp[1];
      r_rel[1] = sin_domega * temp[0] + cos_domega * temp[1];
    }
    else if (n_dim_ == 3) 
      rotate_3d_vector(domega, r_rel, unit_torque);
    normalize_vector(r_rel,n_dim_);
    for (int i=0; i<n_dim_; ++i) {
      it->orientation_[i] = r_rel[i];
      it->position_[i] = position_[i] + anchor_distance_*it->orientation_[i];
    }
  }
}

void Spindle::Integrate() {
  Translate();
  Rotate();
}

std::vector<Object*> Spindle::GetInteractors() {
  interactors_.clear();
  interactors_.push_back(this);
  for (auto it=filaments_.begin(); it!=filaments_.end(); ++it) {
    auto ix_vec = it->GetInteractors();
    interactors_.insert(interactors_.end(), ix_vec.begin(), ix_vec.end());
  }
  return interactors_;
}

void Spindle::Draw(std::vector<graph_struct*> * graph_array) {
  Object::Draw(graph_array);
  for (auto fil=filaments_.begin(); fil!= filaments_.end(); ++fil) {
    fil->Draw(graph_array);
  }
}

