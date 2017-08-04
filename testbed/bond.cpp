#include "bond.h"

/**************************
** Bond member functions **
**************************/
void Bond::Init(Site * s1, Site * s2) {
  s1->AddBond(this,OUTGOING);
  s2->AddBond(this,INCOMING);
  sites_[0] = s1;
  sites_[1] = s2;
  double const * const r1 = s1->GetPosition();
  double const * const r2 = s2->GetPosition();
  diameter_ = s1->GetDiameter();
  length_ = 0;
  for (int i=0;i<n_dim_;++i) {
    orientation_[i] = r2[i] - r1[i];
    length_ += orientation_[i]*orientation_[i];
  }
  length_ = sqrt(length_);
  equil_length_ = length_;
  for (int i=0;i<n_dim_;++i) {
    position_[i] = r1[i] + 0.5*orientation_[i];
    orientation_[i]/=length_;
  }
}

void Bond::ReInit() {
  double const * const r1 = sites_[0]->GetPosition();
  double const * const r2 = sites_[1]->GetPosition();
  diameter_ = sites_[0]->GetDiameter();
  length_ = 0;
  for (int i=0;i<n_dim_;++i) {
    orientation_[i] = r2[i] - r1[i];
    length_ += orientation_[i]*orientation_[i];
  }
  length_ = sqrt(length_);
  for (int i=0;i<n_dim_;++i) {
    position_[i] = r1[i] + 0.5*orientation_[i];
    orientation_[i]/=length_;
  }
}

Site * Bond::GetSite(int i) {
  if (i<0 || i>1) {
    std::cerr << "ERROR! Requested adjacent site out of bounds!\n";
  }
  return sites_[i];
}
Bond * Bond::GetNeighborBond(int i) {
  if (i<0 || i>1) {
    std::cerr << "ERROR! Requested neighboring bond out of bounds!\n";
  }
  return sites_[i]->GetOtherBond(GetOID());
}

directed_bond Bond::GetNeighborDirectedBond(int i) {
  if (i<0 || i>1) {
    std::cerr << "ERROR! Requested neighboring bond out of bounds!\n";
  }
  return sites_[i]->GetOtherDirectedBond(GetOID());
}

void Bond::Report() {
  fprintf(stderr, "  Bond:\n");
  Object::Report();
}
void Bond::ReportSites() {
  Report();
  fprintf(stderr,"    Reporting sites:\n");
  for (int i=0;i<2;++i) {
    sites_[i]->Report();
  }
}

//void Bond::UpdatePosition() {
  //ZeroForce();
  //Integrate();
  ////UpdateSitePositions();
  //// Update end site positions for tracking trajectory for neighbors
  ////UpdateBondPositions();
  ////for (auto bond=v_elements_.begin(); bond!= v_elements_.end(); ++bond) 
    ////bond->UpdatePeriodic();
//}

//[> Integration scheme taken from Yu-Guo Tao,
   //J Chem Phys 122 244903 (2005)
   //Explicit calculation of the friction tensor acting on force vector,

   //r(t+dt) = r(t) + (Xi^-1 . F_s(t)) * dt + dr(t),

   //where friction tensor Xi = friction_par * |u><u| + friction_perp * (I - |u><u|),
   //u is the orientation, F_s is the force from interactions, and dr(t) is the
   //random displacement due to random force,

   //dr(t) = Xi^-1 . F_r * dt,

   //which is treated separately as a random displacement with std dev
   //sqrt(2*kT*dt/gamma_(par/perp)) along par/perp unit vectors
   //relative to rod. */
//void Bond::Integrate() {
  ////Explicit calculation of Xi.F_s
  //for (int i=0; i<n_dim_; ++i) {
    //for (int j=0; j<n_dim_; ++j) {
      //position_[i] += friction_par_*orientation_[i]*orientation_[j]*force_[j]*delta_;
    //}
    //position_[i] += force_[i]*friction_perp_*delta_;
  //}
  ////Add the random displacement dr(t)
  //AddRandomDisplacement();
  ////Update the orientation due to torques and random rotation
  //UpdateOrientation();
//}

//[> Calculates body frame, which returns the vector(s) orthogonal
   //to u(t), then applies random displacements along each
   //orthogonal vector and along u(t) pulled from a distribution
   //with std dev sqrt(2*kT*dt/gamma) where gamma is the friction
   //coefficient along that direction */
//void Bond::AddRandomDisplacement() {
  //// Get vector(s) orthogonal to orientation
  //GetBodyFrame();
  //// First handle the parallel component
  //double mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_par_);
  //for (int i=0; i<n_dim_; ++i)
    //position_[i] += mag * orientation_[i];
  //// Then the perpendicular component(s)
  //for (int j=0; j<n_dim_-1; ++j) {
    //mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_perp_);
    //for (int i=0; i<n_dim_; ++i)
      //position_[i] += mag * body_frame_[n_dim_*j+i];
  //}
  //// Handle the random orientation update after updating orientation from
  //// interaction torques
//}


//[> The orientation update is also from Yu-Guo Tao,

   //u(t+dt) = u(t) + friction_rot^-1 * T_s(t) x u(t) * dt + du(t)

   //where similar to above, du(t) is the reorientation due to
   //random forces, and is treated as random displacement vector(s)
   //orthogonal to u(t) with std dev sqrt(2*kT*dt/friction_rot) */
//void Bond::UpdateOrientation() {
  //// First handle reorientation due to external torques
  //double du[3];
  //cross_product(torque_, orientation_, du, 3); // ndim=3 since torques
  //for (int i=0; i<n_dim_; ++i)
    //orientation_[i] += du[i]*delta_/friction_rot_;
  //// Now handle the random orientation update
    //for (int j=0; j<n_dim_-1; ++j) {
      //double mag = gsl_ran_gaussian_ziggurat(rng_.r, rand_sigma_rot_);
      //for (int i=0; i<n_dim_; ++i)
        //orientation_[i] += mag * body_frame_[n_dim_*j+i];
    //}
  //normalize_vector(orientation_, n_dim_);
//}

//[> calculates vector(s) orthogonal to orientation of rod <]
//void Bond::GetBodyFrame() {
  //if (n_dim_==2) {
    //body_frame_[0] = orientation_[1];
    //body_frame_[1] = -orientation_[0];
  //}
  //else {
    //double vect1[3] = {1.0, 0.0, 0.0};
    //double vect2[3] = {0.0, 1.0, 0.0};
    //if (1.0 - ABS(orientation_[0]) > 1e-2)
      //cross_product(orientation_, vect1, &(body_frame_[0]), n_dim_);
    //else
      //cross_product(orientation_, vect2, &(body_frame_[0]), n_dim_);
    //normalize_vector(&(body_frame_[0]),n_dim_);
    //cross_product(orientation_, &(body_frame_[0]), &(body_frame_[3]), n_dim_);
  //}
//}


