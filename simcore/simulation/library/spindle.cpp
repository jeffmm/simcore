#include "spindle.hpp"

Spindle::Spindle() : Spherocylinder() {}

void Spindle::SetParameters() {
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  diameter_ = sparams_->diameter;
  length_ = sparams_->length;
  n_filaments_bud_ = sparams_->n_filaments_bud;
  n_filaments_mother_ = sparams_->n_filaments_mother;
  k_spring_ = sparams_->k_spring;
  k_align_ = sparams_->k_align;
  spring_length_ = sparams_->spring_length;
  alignment_potential_ = (sparams_->alignment_potential ? true : false);
  anchor_distance_ =
      0.5 * diameter_ + (0.5 + 1e-3) * 1; // FIXME 1 IS FILAMENT DIAMETER
  spb_diameter_ = sparams_->spb_diameter;
  // Force spherocylinder insertion options to match spindle insertion options
  // FIXME
  is_midstep_ = true; // Always true since we have filaments
  SetDiffusion();
}

void Spindle::Init(spindle_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
  filaments_.reserve(n_filaments_bud_ + n_filaments_mother_);
  anchors_.resize(n_filaments_bud_ + n_filaments_mother_);
  bool out_of_bounds;
  int n_insert = 0;
  InsertSpherocylinder();
  GetBodyFrame();
  GenerateAnchorSites();
  filaments_.clear();
  for (int i = 0; i < n_filaments_bud_; ++i) {
    InsertFilament(i);
  }
  for (int i = n_filaments_bud_; i < n_filaments_mother_ + n_filaments_bud_;
       ++i) {
    InsertFilament(i);
  }
}

// Returns true if successful, false otherwise
bool Spindle::InsertFilament(int i) {
  int n_insert = 0;
  bool out_of_bounds;
  Filament fil;
  filaments_.push_back(fil);
  return true;
}

int Spindle::GetCount() {
  int count = 1;
  for (filament_iterator it = filaments_.begin(); it != filaments_.end();
       ++it) {
    count += it->GetCount();
  }
  return count;
}

void Spindle::ZeroForce() {
  Spherocylinder::ZeroForce();
  for (filament_iterator it = filaments_.begin(); it != filaments_.end();
       ++it) {
    it->ZeroForce();
  }
  // FIXME
  // for (anchor_iterator it = anchors_.begin(); it != anchors_.end(); ++it) {
  // it->ZeroForce();
  //}
}

void Spindle::GenerateAnchorSites() {
  // FIXME
  for (int i_fil = 0; i_fil < n_filaments_bud_ + n_filaments_mother_; ++i_fil) {
    // anchors_[i_fil].theta_ = atan(
    // 2.0 * spb_diameter_ * (gsl_rng_uniform_pos(rng_.r) - 0.5) / diameter_);
    // anchors_[i_fil].phi_ = atan(
    // 2.0 * spb_diameter_ * (gsl_rng_uniform_pos(rng_.r) - 0.5) / diameter_);
    // anchors_[i_fil].k_spring_ = k_spring_;
    // anchors_[i_fil].k_align_ = k_align_;
    // anchors_[i_fil].spring_length_ = spring_length_;
    // anchors_[i_fil].alignment_potential_ = alignment_potential_;
  }
  ResetAnchorPositions();
}

void Spindle::ResetAnchorPositions() {
  for (int i_fil = 0; i_fil < n_filaments_bud_ + n_filaments_mother_; ++i_fil) {
    // FIXME
    // std::copy(orientation_, orientation_ + 3, anchors_[i_fil].orientation_);
    if (i_fil >= n_filaments_bud_) {
      for (int i = 0; i < n_dim_; ++i) {
        // FIXME
        // anchors_[i_fil].orientation_[i] *= -1;
      }
    }
    // FIXME
    // rotate_vector(anchors_[i_fil].orientation_, &body_frame_[0],
    // anchors_[i_fil].theta_);
    // rotate_vector(anchors_[i_fil].orientation_, &body_frame_[3],
    // anchors_[i_fil].phi_);
    int sign = (i_fil >= n_filaments_bud_ ? -1 : 1);
    for (int i = 0; i < n_dim_; ++i) {
      // FIXME
      // anchors_[i_fil].position_[i] =
      // position_[i] + sign * 0.5 * orientation_[i] * length_ +
      // FIXME
      // anchor_distance_ * anchors_[i_fil].orientation_[i];
    }
  }
}

void Spindle::UpdatePosition(bool midstep) {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<std::vector<Filament>::iterator,
                        std::vector<Filament>::iterator>>
      chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = filaments_.size() / max_threads;
  filament_iterator cur_iter = filaments_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    filament_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, filaments_.end()));
#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto it = chunks[i].first; it != chunks[i].second; ++it) {
        it->UpdatePosition(midstep);
      }
    }
  }
#else
  for (filament_iterator it = filaments_.begin(); it != filaments_.end();
       ++it) {
    it->UpdatePosition(midstep);
  }
#endif
  SetPrevPosition(position_);
  ApplyForcesTorques();
  Integrate();
  UpdatePeriodic();
}

void Spindle::ApplyForcesTorques() {
  // FIXME
  //for (anchor_iterator it = anchors_.begin(); it != anchors_.end(); ++it) {
    //// First calculate translational forces
    //// double f_mag = dot_product(n_dim_, it->orientation_, it->force_);
    //double dr[3] = {0, 0, 0};
    //for (int i = 0; i < n_dim_; ++i) {
      //// force_[i] += f_mag * it->orientation_[i];
      //force_[i] += it->force_[i];
      //dr[i] = it->position_[i] - position_[i];
    //}
    //// Then calculate torques
    //double i_torque[3] = {0, 0, 0};
    //cross_product(dr, it->force_, i_torque, 3);
    //// 3-dimensions due to torques
    //for (int i = 0; i < 3; ++i) {
      //torque_[i] += i_torque[i];
      //// And add torque from alignment potential
      //if (alignment_potential_) {
        //torque_[i] += it->torque_[i];
      //}
    //}
  //}
}

void Spindle::Integrate() {
  Spherocylinder::Integrate();
  // Now I need to update the positions of the anchors based on the new
  // positions and orientations using their euler angles
  ResetAnchorPositions();
}

void Spindle::GetInteractors(std::vector<Object *> &ix) {
  ix.push_back(this);
  for (auto it = filaments_.begin(); it != filaments_.end(); ++it) {
    it->GetInteractors(ix);
  }
}

void Spindle::Draw(std::vector<graph_struct *> &graph_array) {
  Object::Draw(graph_array);
  for (auto fil = filaments_.begin(); fil != filaments_.end(); ++fil) {
    fil->Draw(graph_array);
  }
}
