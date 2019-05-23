#include "spindle.h"

Spindle::Spindle() : Spherocylinder() {
  color_ = params_->spindle.color;
  draw_ = draw_type::_from_string(params_->spindle.draw_type.c_str());
  diameter_ = params_->spindle.diameter;
  length_ = params_->spindle.length;
  n_filaments_bud_ = params_->spindle.n_filaments_bud;
  n_filaments_mother_ = params_->spindle.n_filaments_mother;
  k_spring_ = params_->spindle.k_spring;
  k_align_ = params_->spindle.k_align;
  spring_length_ = params_->spindle.spring_length;
  alignment_potential_ = (params_->spindle.alignment_potential ? true : false);
  anchor_distance_ = 0.5*diameter_ + (0.5 + 1e-3)*params_->filament.diameter;
  spb_diameter_ = params_->spindle.spb_diameter;
  // Force spherocylinder insertion options to match spindle insertion options
  params_->spherocylinder.insertion_type = params_->spindle.insertion_type;
  is_midstep_ = true; // Always true since we have filaments
  SetDiffusion();
}

void Spindle::Init() {
  filaments_.reserve(n_filaments_bud_+n_filaments_mother_);
  anchors_.resize(n_filaments_bud_+n_filaments_mother_);
  bool out_of_bounds;
  int n_insert = 0;
  do {
    InsertSpherocylinder();
    GetBodyFrame();
    double buffer = MIN(1.5*params_->filament.min_length,1.1*params_->filament.max_length);
    if (out_of_bounds = CheckBounds(buffer)) continue;
    GenerateAnchorSites();
    filaments_.clear();
    for (int i=0;i<n_filaments_bud_;++i) {
      if (out_of_bounds = !InsertFilament(i)) {
        break;
      }
    }
    if (out_of_bounds) continue;
    for (int i=n_filaments_bud_;i<n_filaments_mother_+n_filaments_bud_;++i) {
      if (out_of_bounds = !InsertFilament(i)) {
        break;
      }
    }
  } while (out_of_bounds);
}

// Returns true if successful, false otherwise
bool Spindle::InsertFilament(int i) {
  int n_insert = 0;
  bool out_of_bounds;
  do {
    Filament fil;
    filaments_.push_back(fil);
    filaments_.back().SetAnchor(&anchors_[i]);
    filaments_.back().Init(true);
    if (out_of_bounds = filaments_.back().CheckBounds()) {
      filaments_.pop_back();
    }
    if (out_of_bounds && n_insert++ > 1000) {
      break;
    }
  } while (out_of_bounds);
  return !out_of_bounds;
}

int Spindle::GetCount() {
  int count=1;
  for (filament_iterator it=filaments_.begin();it!=filaments_.end();++it) {
    count += it->GetCount();
  }
  return count;
}

void Spindle::ZeroForce() {
  Spherocylinder::ZeroForce();
  for (filament_iterator it=filaments_.begin();it!=filaments_.end();++it) {
    it->ZeroForce();
  }
  for (anchor_iterator it=anchors_.begin(); it!=anchors_.end(); ++it) {
    it->ZeroForce();
  }
}

void Spindle::GenerateAnchorSites() {
  for (int i_fil=0;i_fil<n_filaments_bud_+n_filaments_mother_;++i_fil) {
    anchors_[i_fil].theta_ = atan(2.0*spb_diameter_*(gsl_rng_uniform_pos(rng_.r)-0.5)/diameter_);
    anchors_[i_fil].phi_ = atan(2.0*spb_diameter_*(gsl_rng_uniform_pos(rng_.r)-0.5)/diameter_);
    anchors_[i_fil].k_spring_ = k_spring_;
    anchors_[i_fil].k_align_ = k_align_;
    anchors_[i_fil].spring_length_ = spring_length_;
    anchors_[i_fil].alignment_potential_ = alignment_potential_;
  }
  ResetAnchorPositions();
}

void Spindle::ResetAnchorPositions() {
  for (int i_fil=0;i_fil<n_filaments_bud_+n_filaments_mother_;++i_fil) {
    std::copy(orientation_,orientation_+3,anchors_[i_fil].orientation_);
    if (i_fil>=n_filaments_bud_) {
      for (int i=0;i<n_dim_;++i) {
        anchors_[i_fil].orientation_[i] *= -1;
      }
    }
    rotate_vector(anchors_[i_fil].orientation_, &body_frame_[0], anchors_[i_fil].theta_);
    rotate_vector(anchors_[i_fil].orientation_, &body_frame_[3], anchors_[i_fil].phi_);
    int sign = (i_fil >= n_filaments_bud_ ? -1 : 1);
    for (int i=0;i<n_dim_;++i) {
      anchors_[i_fil].position_[i] = position_[i] + sign*0.5*orientation_[i]*length_ + anchor_distance_*anchors_[i_fil].orientation_[i];
    }
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
  SetPrevPosition(position_);
  ApplyForcesTorques();
  Integrate();
  UpdatePeriodic();
}

void Spindle::ApplyForcesTorques() {
  for (anchor_iterator it=anchors_.begin(); it!=anchors_.end(); ++it) {
    // First calculate translational forces
    //double f_mag = dot_product(n_dim_, it->orientation_, it->force_);
    double dr[3] = {0,0,0};
    for (int i=0; i<n_dim_; ++i) {
      //force_[i] += f_mag * it->orientation_[i];
      force_[i] += it->force_[i];
      dr[i] = it->position_[i]-position_[i];
    }
    // Then calculate torques
    double i_torque[3] = {0,0,0};
    cross_product(dr, it->force_, i_torque, 3);
     //3-dimensions due to torques
    for (int i=0; i<3; ++i) {
      torque_[i] += i_torque[i];
       //And add torque from alignment potential
      if (alignment_potential_) {
        torque_[i] += it->torque_[i];
      }
    }
  }
}

void Spindle::Integrate() {
  Spherocylinder::Integrate();
  // Now I need to update the positions of the anchors based on the new
  // positions and orientations using their euler angles
  ResetAnchorPositions();
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

