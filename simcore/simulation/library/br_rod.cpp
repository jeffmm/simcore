#include "br_rod.hpp"

BrRod::BrRod(unsigned long seed) : Object(seed) {
  std::fill(body_frame_, body_frame_ + 6, 0.0);
}

void BrRod::InsertRod(std::string insertion_type, double buffer) {
  if (buffer < 0) {
    buffer = length_;
  }
  if (insertion_type.compare("random") == 0) {
    InsertRandom();
  } else if (insertion_type.compare("random_oriented") == 0) {
    InsertRandom();
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = 1.0;
  } else if (insertion_type.compare("centered_random") == 0) {
    std::fill(position_, position_ + 3, 0.0);
    rng_.RandomUnitVector(n_dim_, orientation_);
  } else if (insertion_type.compare("centered_oriented") == 0) {
    std::fill(position_, position_ + 3, 0.0);
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = 1.0;
  } else {
    Logger::Error("BrRod insertion type not recognized!");
  }
}

/* Integration scheme taken from Yu-Guo Tao,
   J Chem Phys 122 244903 (2005)
   Explicit calculation of the friction tensor acting on force vector,

   r(t+dt) = r(t) + (Xi^-1 . F_s(t)) * dt + dr(t),

   where friction tensor Xi = gamma_par * |u><u| + gamma_perp * (I - |u><u|),
   u is the orientation, F_s is the force from interactions, and dr(t) is the
   random displacement due to random force,

   dr(t) = Xi^-1 . F_r * dt,

   which is treated separately as a random displacement with std dev
   sqrt(2*kT*dt/gamma_(par/perp)) along par/perp unit vectors
   relative to rod. */
void BrRod::Integrate() {
  // Explicit calculation of Xi.F_s
  double dr[3] = {0, 0, 0};
  for (int i = 0; i < n_dim_; ++i) {
    for (int j = 0; j < n_dim_; ++j) {
      dr[i] +=
          gamma_par_ * orientation_[i] * orientation_[j] * force_[j] * delta_;
    }
    dr[i] += force_[i] * gamma_perp_ * delta_;
    position_[i] += dr[i];
  }
  // Reorientation due to external torques
  double du[3];
  cross_product(torque_, orientation_, du, 3); // ndim=3 since torques
  for (int i = 0; i < n_dim_; ++i) {
    orientation_[i] += du[i] * delta_ / gamma_rot_;
  }
  // Add the random displacement dr(t)
  AddRandomDisplacement();
  // Update the orientation due to torques and random rotation
  AddRandomReorientation();
  UpdatePeriodic();
}

/* Calculates body frame, which returns the vector(s) orthogonal
   to u(t), then applies random displacements along each
   orthogonal vector and along u(t) pulled from a distribution
   with std dev sqrt(2*kT*dt/gamma) where gamma is the friction
   coefficient along that direction */
void BrRod::AddRandomDisplacement() {
  // Get vector(s) orthogonal to orientation
  GetBodyFrame();
  // First handle the parallel component
  double mag = rng_.RandomNormal(diffusion_par_);
  for (int i = 0; i < n_dim_; ++i)
    position_[i] += mag * orientation_[i];
  // Then the perpendicular component(s)
  for (int j = 0; j < n_dim_ - 1; ++j) {
    mag = rng_.RandomNormal(diffusion_perp_);
    for (int i = 0; i < n_dim_; ++i)
      position_[i] += mag * body_frame_[n_dim_ * j + i];
  }
  // Handle the random orientation update after updating orientation from
  // interaction torques
}

/* The orientation update is also from Yu-Guo Tao,

   u(t+dt) = u(t) + gamma_rot^-1 * T_s(t) x u(t) * dt + du(t)

   where similar to above, du(t) is the reorientation due to
   random forces, and is treated as random displacement vector(s)
   orthogonal to u(t) with std dev sqrt(2*kT*dt/gamma_rot) */
void BrRod::AddRandomReorientation() {
  // Now handle the random orientation update
  for (int j = 0; j < n_dim_ - 1; ++j) {
    double mag = rng_.RandomNormal(diffusion_rot_);
    for (int i = 0; i < n_dim_; ++i) {
      orientation_[i] += mag * body_frame_[n_dim_ * j + i];
    }
  }
  normalize_vector(orientation_, n_dim_);
}

void BrRod::SetDiffusion() {
  // Sets the diffusion in accordance to Lowen, Phys. Rev. E, 1994.
  double L = length_ + diameter_;
  double p = L / diameter_;
  double log_p = log(p);
  gamma_par_ = 2.0 * L / 3.0 / (log_p - 0.207 + 0.980 / p - 0.133 / SQR(p));
  gamma_perp_ = 4.0 * L / 3.0 / (log_p + 0.839 + 0.185 / p + 0.233 / SQR(p));
  gamma_rot_ = CUBE(L) / 9.0 / (log_p - 0.662 + 0.917 / p - 0.050 / SQR(p));
  diffusion_par_ = sqrt(2 * delta_ / gamma_par_);
  diffusion_perp_ = sqrt(2 * delta_ / gamma_perp_);
  diffusion_rot_ = sqrt(2 * delta_ / gamma_rot_);
}

void BrRod::GetBodyFrame() {
  if (n_dim_ == 2) {
    body_frame_[0] = orientation_[1];
    body_frame_[1] = -orientation_[0];
  } else {
    double vect1[3] = {1.0, 0.0, 0.0};
    double vect2[3] = {0.0, 1.0, 0.0};
    if (1.0 - ABS(orientation_[0]) > 1e-2)
      cross_product(orientation_, vect1, &(body_frame_[0]), n_dim_);
    else
      cross_product(orientation_, vect2, &(body_frame_[0]), n_dim_);
    normalize_vector(&(body_frame_[0]), n_dim_);
    cross_product(orientation_, &(body_frame_[0]), &(body_frame_[3]), n_dim_);
  }
}
