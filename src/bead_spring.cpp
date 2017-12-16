#include "bead_spring.h"

BeadSpring::BeadSpring() : Mesh() {
  SetParameters();
} 

void BeadSpring::SetParameters() {
  color_ = params_->bead_spring.color;
  draw_ = draw_type::_from_string(params_->bead_spring.draw_type.c_str());
  length_ = params_->bead_spring.length;
  persistence_length_ = params_->bead_spring.persistence_length;
  diameter_ = params_->bead_spring.diameter;
  max_bond_length_ = params_->bead_spring.max_bond_length;
  bond_spring_ = params_->bead_spring.bond_spring;
  driving_factor_ = params_->bead_spring.driving_factor;
  stoch_flag_ = params_->stoch_flag; // determines whether we are using stochastic forces
  eq_steps_ = params_->n_steps_equil;
  eq_steps_count_ = 0;
}

bool BeadSpring::CheckBounds(double buffer) {
  for (auto site=sites_.begin();site!=sites_.end(); ++site) {
    if (site->CheckBounds()) return true;
  }
  return false;
}

void BeadSpring::Init(bool force_overlap) {
  InitElements();
  InsertBeadSpring(force_overlap);
  UpdatePrevPositions();
  CalculateAngles();
  SetDiffusion();
}

void BeadSpring::InitElements() {
  n_bonds_max_ = (int) ceil(length_/(0.8*max_bond_length_));
  Reserve(n_bonds_max_);
  //Allocate control structures
  int n_sites_max = n_bonds_max_+1;
  cos_thetas_.resize(n_sites_max-2); //max_sites-2
}

void BeadSpring::GenerateProbableOrientation() {
  /* This updates the current orientation with a generated probable orientation where
  we generate random theta pulled from probability distribution P(th) = exp(k cos(th))
  where k is the persistence_length of the filament
  If k is too large, there is enormous imprecision in this calculation since sinh(k)
  is very large so to fix this I introduce an approximate distribution that is valid 
  for large k */
  double theta;
  if (persistence_length_ == 0) {
    theta = gsl_rng_uniform_pos(rng_.r) * M_PI;
  }
  else if (persistence_length_ < 100) {
    theta = gsl_rng_uniform_pos(rng_.r) * M_PI;
    theta = acos( log( exp(-persistence_length_/bond_length_) + 
          2.0*gsl_rng_uniform_pos(rng_.r)*sinh(persistence_length_/bond_length_) ) 
        / (persistence_length_/bond_length_) );
  }
  else {
    theta = acos( (log( 2.0*gsl_rng_uniform_pos(rng_.r)) - 
          log(2.0) + persistence_length_/bond_length_)
        /(persistence_length_/bond_length_) );
  }
  double new_orientation[3] = {0, 0, 0};
  if (n_dim_==2) {
    theta = (gsl_rng_uniform_int(rng_.r,2)==0 ? -1 : 1) * theta;
    new_orientation[0] = cos(theta);
    new_orientation[1] = sin(theta);
  }
  else {
    double phi = gsl_rng_uniform_pos(rng_.r) * 2.0 * M_PI;
    new_orientation[0] = sin(theta)*cos(phi);
    new_orientation[1] = sin(theta)*sin(phi);
    new_orientation[2] = cos(theta);
  }
  rotate_orientation_vector(n_dim_, new_orientation, orientation_);
  std::copy(new_orientation,new_orientation+3,orientation_);
}

void BeadSpring::InsertFirstBond() {
  if (params_->bead_spring.insertion_type.compare("random") == 0) {
    InitRandomSite(diameter_);
    AddRandomBondToTip(bond_length_);
  }
  else if (params_->bead_spring.insertion_type.compare("random_oriented") == 0) {
    InitRandomSite(diameter_);
    std::fill(orientation_,orientation_+3,0.0);
    orientation_[n_dim_-1] = (gsl_rng_uniform_pos(rng_.r) > 0.5 ? 1.0 : -1.0);
    AddBondToTip(orientation_, bond_length_);
  }
  else if (params_->bead_spring.insertion_type.compare("centered_oriented") == 0) {
    std::fill(position_,position_+3,0.0);
    position_[n_dim_-1] = -0.5*length_;
    InitSiteAt(position_,diameter_);
    std::fill(orientation_,orientation_+3,0.0);
    orientation_[n_dim_-1] = 1.0;
    AddBondToTip(orientation_, bond_length_);
  }
  else if (params_->bead_spring.insertion_type.compare("centered_random") == 0) {
    generate_random_unit_vector(n_dim_,orientation_,rng_.r);
    for (int i=0;i<n_dim_;++i) {
      position_[i] = - 0.5*length_*orientation_[i];
    }
    InitSiteAt(position_,diameter_);
    AddBondToTip(orientation_, bond_length_);
  }
}

void BeadSpring::InsertBeadSpring(bool force_overlap) {
  bool out_of_bounds = true;
  do {
    out_of_bounds = false;
    bond_length_ = 0.8*max_bond_length_;
    Clear();
    InsertFirstBond();
    if (!force_overlap && (out_of_bounds = sites_[n_sites_-1].CheckBounds())) continue;
    SetOrientation(bonds_[n_bonds_-1].GetOrientation());
    bool probable_orientation = (params_->bead_spring.insertion_type.compare("simple_crystal") != 0 && params_->bead_spring.insertion_type.compare("random_oriented") != 0);
    for (int i=0;i<n_bonds_max_-1;++i) {
      if (probable_orientation) {
        GenerateProbableOrientation();
      }
      AddBondToTip(orientation_, bond_length_);
      if (!force_overlap && (out_of_bounds = sites_[n_sites_-1].CheckBounds())) break;
    }
  } while (out_of_bounds);
  for (bond_iterator bond=bonds_.begin(); bond!=bonds_.end(); ++bond) {
    bond->SetColor(color_,draw_);
  }
  UpdateBondPositions();
  UpdateSiteOrientations();
}

void BeadSpring::InsertAt(double *pos, double *u) {
  bond_length_ = 0.8*max_bond_length_;
  for (int i=0;i<n_dim_; ++i) {
    position_[i] = pos[i] - 0.5*length_*u[i];
  }
  InitSiteAt(position_,diameter_);
  std::copy(u,u+3,orientation_);
  AddBondToTip(orientation_, bond_length_);
  if (sites_[n_sites_-1].CheckBounds()) {
    error_exit("BeadSpring inserted manually out of bounds.");
  }
  SetOrientation(bonds_[n_bonds_-1].GetOrientation());
  for (int i=0;i<n_bonds_max_-1;++i) {
    if (params_->bead_spring.insertion_type.compare("simple_crystal") != 0) {
      GenerateProbableOrientation();
    }
    AddBondToTip(orientation_, bond_length_);
    if (sites_[n_sites_-1].CheckBounds()) {
      error_exit("BeadSpring inserted manually out of bounds.");
    }
  }
  UpdateBondPositions();
  UpdateSiteOrientations();
  UpdatePrevPositions();
  CalculateAngles();
  SetDiffusion();
}

void BeadSpring::SetDiffusion() {
  rand_sigma_ = sqrt(24.0*diameter_/delta_);
}

double const BeadSpring::GetVolume() {
  if (n_dim_ == 2) {
    return diameter_*length_ + 0.25*M_PI*diameter_*diameter_;
  }
  if (n_dim_ == 3) {
    return 0.25*M_PI*diameter_*diameter_*length_ + 1.0/6.0 * M_PI * diameter_*diameter_*diameter_;
  }
}

void BeadSpring::UpdatePosition() {
  //ApplyForcesTorques();
  Integrate();
  eq_steps_count_++;
}

/*******************************************************************************
  BD algorithm for inextensible wormlike chains with anisotropic friction
  Montesi, Morse, Pasquali. J Chem Phys 122, 084903 (2005).
********************************************************************************/
void BeadSpring::Integrate() {
  CalculateAngles();
  CalculateTangents();
  AddRandomForces();
  AddBendingForces();
  AddSpringForces();
  UpdateSitePositions();
  UpdateBondPositions();
  UpdateSiteOrientations();
}

void BeadSpring::UpdateSiteOrientations() {
  for (int i=0; i<n_sites_-1; ++i) {
    sites_[i].SetOrientation(bonds_[i].GetOrientation());
  }
  sites_[n_sites_-1].SetOrientation(bonds_[n_bonds_-1].GetOrientation());
}

void BeadSpring::CalculateAngles() {
  double cos_angle;
  for (int i_site=0; i_site<n_sites_-2; ++i_site) {
    double const * const u1 = sites_[i_site].GetOrientation();
    double const * const u2 = sites_[i_site+1].GetOrientation();
    cos_angle = dot_product(n_dim_, u1, u2);
    cos_thetas_[i_site] = cos_angle;
  }
}

void BeadSpring::CalculateTangents() {
  for (auto it = sites_.begin(); it!=sites_.end(); ++it) {
    it->CalcTangent();
  }
}

void BeadSpring::AddRandomForces() {
  if (!stoch_flag_) return;
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    for (int i=0; i<n_dim_; ++i) {
      double kick = gsl_rng_uniform_pos(rng_.r) - 0.5;
      force_[i] = kick*rand_sigma_;
    }
    sites_[i_site].AddForce(force_);
  }
}

void BeadSpring::AddBendingForces() {
  /* These forces are calculated from the potential U = k(theta-pi)^2 where
   * theta is the angle between two neighboring bonds consisting of three
   * beads.  I will include the derivation of the site forces in the
   * documentation (JMM 20171215) */
  for (int i_theta=0; i_theta<n_bonds_-1; ++i_theta) {
    double udotu = dot_product(n_dim_, sites_[i_theta].GetOrientation(), 
        sites_[i_theta+1].GetOrientation());
    double kappa;
    if (udotu - 1 < 1e-10) {
      kappa = 0;
    }
    else {
      kappa = 2*persistence_length_*(acos(cos_thetas_[i_theta])-
        M_PI)/sqrt(1-SQR(cos_thetas_[i_theta]));
    }
    double rho1 = 0;
    double rho2 = 0;
    for (int i=0;i<n_dim_;++i) {
      double temp = sites_[i_theta+1].GetPosition()[i] - sites_[i_theta].GetPosition()[i];
      rho1 += temp*temp;
      temp = sites_[i_theta+2].GetPosition()[i] - sites_[i_theta+1].GetPosition()[i];
      rho2 += temp*temp;
    }
    rho1 = 1.0/sqrt(rho1);
    rho2 = 1.0/sqrt(rho2);
    double const * const u1 = sites_[i_theta].GetOrientation();
    double const * const u2 = sites_[i_theta+1].GetOrientation();
    // Forces for the first of three sites, which is part of the force on the second site
    for (int i=0; i<n_dim_; ++i) {
      force_[i] = kappa*(-u2[i]*rho1 + udotu*u1[i]*rho1);
    }
    sites_[i_theta].AddForce(force_);
    sites_[i_theta+1].SubForce(force_);
    // Forces for the third of three sites, which is part of the force on the second site
    for (int i=0; i<n_dim_; ++i) {
      force_[i] = kappa*(u1[i]*rho2 - udotu*u2[i]*rho2);
    }
    sites_[i_theta+2].AddForce(force_);
    sites_[i_theta+1].SubForce(force_);
  }
}

void BeadSpring::AddSpringForces() {
  /* FENE springs have the potential U(r) = -1/2 k r0 ln(1 - (r/r0)^2) */
  for(int i_bond=0; i_bond<n_bonds_; ++i_bond) {
    double r = bonds_[i_bond].GetLength();
    double f_mag = bond_spring_ * (r-max_bond_length_);
    for (int i=0;i<n_dim_;++i) {
      force_[i] = f_mag * sites_[i_bond].GetOrientation()[i];
    }
    sites_[i_bond+1].SubForce(force_);
    sites_[i_bond].AddForce(force_);
  }

  //for(int i_bond=0; i_bond<n_bonds_; ++i_bond) {
    //double r = bonds_[i_bond].GetLength();
    //double f_mag = bond_spring_ * r / (1 - SQR(r)/SQR(max_bond_length_));
    //for (int i=0;i<n_dim_;++i) {
      //force_[i] = f_mag * sites_[i_bond].GetOrientation()[i];
    //}
    //sites_[i_bond+1].SubForce(force_);
    //sites_[i_bond].AddForce(force_);
  //}
}

void BeadSpring::UpdateSitePositions() {
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    double const * const f_site = sites_[i_site].GetForce();
    double const * const r_site = sites_[i_site].GetPosition();
    sites_[i_site].SetPrevPosition(r_site);
    for (int i=0; i<n_dim_; ++i) {
      position_[i] = r_site[i] + f_site[i] * delta_ / diameter_;
    }
    sites_[i_site].SetPosition(position_);
  }
  // Next, update orientation vectors
  double u_mag, r_diff[3];
  for (int i_site=0; i_site<n_sites_-1; ++i_site) {
    double const * const r_site1 = sites_[i_site].GetPosition();
    double const * const r_site2 = sites_[i_site+1].GetPosition();
    u_mag = 0.0;
    for (int i=0; i<n_dim_; ++i) {
      r_diff[i] = r_site2[i] - r_site1[i];
      u_mag += SQR(r_diff[i]);
    }
    u_mag = sqrt(u_mag);
    for (int i=0; i<n_dim_; ++i)
      r_diff[i]/=u_mag;
    sites_[i_site].SetOrientation(r_diff);
    sites_[i_site].UpdatePeriodic();
  }
  sites_[n_sites_-1].SetOrientation(sites_[n_sites_-2].GetOrientation());
}

void BeadSpring::UpdateAvgPosition() {
  std::fill(position_, position_+3, 0.0);
  std::fill(orientation_, orientation_+3, 0.0);
  for (auto site_it=sites_.begin(); site_it!=sites_.end(); ++site_it) {
    double const * const site_pos = site_it->GetPosition();
    double const * const site_u = site_it->GetOrientation();
    for (int i=0; i<n_dim_; ++i) {
      position_[i] += site_pos[i];
      orientation_[i] += site_u[i];
    }
  }
  normalize_vector(orientation_, n_dim_);
  for (int i=0; i<n_dim_; ++i) {
    position_[i] /= n_sites_;
  }
}

void BeadSpring::ApplyForcesTorques() {
  ApplyInteractionForces();
}

//FIXME
void BeadSpring::ApplyInteractionForces() {}

void BeadSpring::Draw(std::vector<graph_struct*> * graph_array) {
  for (auto site=sites_.begin(); site!= sites_.end(); ++site) {
    site->Draw(graph_array);
  }
}

// Scale bond and site positions from new unit cell
void BeadSpring::ScalePosition() {
  // scale first bond position using new unit cell
  bonds_[0].ScalePosition();
  // then reposition sites based on first bond position
  // handle first site
  double r[3];
  double const * const bond_r = bonds_[0].GetPosition();
  double const * const bond_u = bonds_[0].GetOrientation();
  for (int i=0; i<n_dim_; ++i)
    r[i] = bond_r[i] - 0.5*bond_length_*bond_u[i];
  sites_[0].SetPosition(r);
  // then handle remaining sites
  for (int i_site=1; i_site<n_sites_; ++i_site) {
    double const * const prev_r = sites_[i_site-1].GetPosition();
    double const * const prev_u = sites_[i_site-1].GetOrientation();
    for (int i=0; i<n_dim_; ++i)
      r[i] = prev_r[i] + bond_length_*prev_u[i];
    sites_[i_site].SetPosition(r);
  }
  // update remaining bond positions
  UpdateBondPositions();
}


void BeadSpring::GetAvgOrientation(double * au) {
  double avg_u[3] = {0.0, 0.0, 0.0};
  int size=0;
  for (auto it=sites_.begin(); it!=sites_.end(); ++it) {
    double const * const u = it->GetOrientation();
    for (int i=0; i<n_dim_; ++i)
      avg_u[i] += u[i];
    size++;
  }
  if (size == 0)
    error_exit("Something went wrong in GetAvgOrientation!");
  for (int i=0; i<n_dim_; ++i)
    avg_u[i]/=size;
  std::copy(avg_u, avg_u+3, au);
}

void BeadSpring::GetAvgPosition(double * ap) {
  double avg_p[3] = {0.0, 0.0, 0.0};
  int size=0;
  for (auto it=sites_.begin(); it!=sites_.end(); ++it) {
    double const * const p = it->GetPosition();
    for (int i=0; i<n_dim_; ++i)
      avg_p[i] += p[i];
    size++;
  }
  if (size == 0)
    error_exit("Something went wrong in GetAvgPosition!");
  for (int i=0; i<n_dim_; ++i)
    avg_p[i]/=size;
  std::copy(avg_p, avg_p+3, ap);
}

void BeadSpring::ReportAll() {
  printf("cos_thetas:\n  {");
  for (int i=0; i<n_sites_-2; ++i)
    printf(" %5.5f ",cos_thetas_[i]);
  printf("}\n");
}

/* The spec output for one bead_spring is:
    diameter
    length
    persistence_length
    n_sites,
    all site positions
    */

void BeadSpring::WriteSpec(std::fstream &ospec){
  ospec.write(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  ospec.write(reinterpret_cast<char*>(&length_), sizeof(length_));
  ospec.write(reinterpret_cast<char*>(&persistence_length_), sizeof(persistence_length_));
  ospec.write(reinterpret_cast<char*>(&n_sites_), sizeof(n_sites_));
  double temp[3];
  for (int i=0; i<n_sites_; ++i) {
    double const * const pos = sites_[i].GetPosition();
    std::copy(pos, pos+3, temp);
    for (auto& p : temp) 
      ospec.write(reinterpret_cast<char*>(&p), sizeof(p));
  }
  return;
}

void BeadSpring::ReadSpec(std::fstream &ispec) {
  if (ispec.eof()) return;
  double r_site[3];
  ispec.read(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  ispec.read(reinterpret_cast<char*>(&length_), sizeof(length_));
  ispec.read(reinterpret_cast<char*>(&persistence_length_), sizeof(persistence_length_));
  ispec.read(reinterpret_cast<char*>(&n_sites_), sizeof(n_sites_));
  sites_.resize(n_sites_, sites_[0]);
  // Get initial site position
  // Initialize sites from relative positions
  for (int i_site=0; i_site<n_sites_; ++i_site) {
    for (int i=0; i<3; ++i)
      ispec.read(reinterpret_cast<char*>(&r_site[i]), sizeof(double));
    sites_[i_site].SetPosition(r_site);
  }
  UpdateBondPositions();
  UpdateSiteOrientations();
  CalculateAngles();
}

void BeadSpring::WritePosit(std::fstream &oposit) {
  double avg_pos[3], avg_u[3];
  GetAvgPosition(avg_pos);
  GetAvgOrientation(avg_u);
  std::copy(avg_pos,avg_pos+3,position_);
  UpdatePeriodic();
  for (auto& pos : position_)
    oposit.write(reinterpret_cast<char*>(&pos), sizeof(pos));
  for (auto& spos : scaled_position_)
    oposit.write(reinterpret_cast<char*>(&spos), sizeof(spos));
  for (auto& u : avg_u) 
    oposit.write(reinterpret_cast<char*>(&u), sizeof(u));
  oposit.write(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  oposit.write(reinterpret_cast<char*>(&length_), sizeof(length_));
}

void BeadSpring::ReadPosit(std::fstream &iposit) {
  if (iposit.eof()) return;
  double avg_pos[3], avg_u[3], s_pos[3];
  for (int i=0; i<3; ++i)
    iposit.read(reinterpret_cast<char*>(&avg_pos[i]), sizeof(double));
  for (int i=0; i<3; ++i)
    iposit.read(reinterpret_cast<char*>(&s_pos[i]), sizeof(double));
  for (int i=0; i<3; ++i)
    iposit.read(reinterpret_cast<char*>(&avg_u[i]), sizeof(double));
  iposit.read(reinterpret_cast<char*>(&diameter_), sizeof(diameter_));
  iposit.read(reinterpret_cast<char*>(&length_), sizeof(length_));
  // Initialize first bond position
  for (int i=0; i<n_dim_; ++i) 
    avg_pos[i] = avg_pos[i] - 0.5*(length_ - bond_length_)*avg_u[i];
  for (int i_bond=0; i_bond<n_bonds_; ++i_bond) {
    bonds_[i_bond].SetPosition(avg_pos);
    bonds_[i_bond].SetOrientation(avg_u);
    bonds_[i_bond].SetDiameter(diameter_);
    bonds_[i_bond].UpdatePeriodic();
    // Set next bond position
    for (int i=0; i<n_dim_; ++i) 
      avg_pos[i] += bond_length_*avg_u[i];
  }
}

void BeadSpring::WriteCheckpoint(std::fstream &ocheck) {
  void * rng_state = gsl_rng_state(rng_.r);
  size_t rng_size = gsl_rng_size(rng_.r);
  ocheck.write(reinterpret_cast<char*>(&rng_size), sizeof(size_t));
  ocheck.write(reinterpret_cast<char*>(rng_state), rng_size);
  WriteSpec(ocheck);
}

void BeadSpring::ReadCheckpoint(std::fstream &icheck) {
  if (icheck.eof()) return;
  void * rng_state = gsl_rng_state(rng_.r);
  size_t rng_size;
  icheck.read(reinterpret_cast<char*>(&rng_size), sizeof(size_t));
  icheck.read(reinterpret_cast<char*>(rng_state), rng_size);
  ReadSpec(icheck);
  n_sites_ = n_bonds_+1;
  double r[3];
  for (int i=0; i<3; ++i)
    r[i] = bonds_[0].GetPosition()[i] - 0.5*bond_length_*bonds_[0].GetOrientation()[i];
  sites_[0].SetPosition(r);
  sites_[0].SetOrientation(bonds_[0].GetOrientation());
  for (int i_bond=1; i_bond<n_bonds_; ++i_bond) {
    for (int i=0; i<3; ++i)
      r[i] += bond_length_*sites_[i_bond-1].GetOrientation()[i];
    sites_[i_bond].SetPosition(r);
    sites_[i_bond].SetOrientation(bonds_[i_bond].GetOrientation());
  }
  for (int i=0; i<3; ++i)
    r[i] += bond_length_*sites_[n_sites_-2].GetOrientation()[i];
  sites_[n_sites_-1].SetPosition(r); 
  sites_[n_sites_-1].SetOrientation(bonds_[n_bonds_-1].GetOrientation());
  UpdateBondPositions();
  //Reallocate control structures
  cos_thetas_.resize(n_sites_-2); //max_sites-2
}

void BeadSpringSpecies::InitAnalysis() {
  time_ = 0;
  if (params_->bead_spring.theta_analysis) {
    InitThetaAnalysis();
  }
  if (params_->bead_spring.lp_analysis) {
    InitMse2eAnalysis();
  }
  RunAnalysis();
}

void BeadSpringSpecies::InitMse2eAnalysis() {
  std::string fname = params_->run_name;
  fname.append("_bead_spring.mse2e");
  mse2e_file_.open(fname, std::ios::out);
  mse2e_file_ << "mse2e_analysis_file\n";
  mse2e_file_ << "length diameter bond_length persistence_length driving ndim nsteps nspec delta theory\n";
  auto it=members_.begin();
  double l = it->GetLength();
  double d = it->GetDiameter();
  double cl = it->GetBondLength();
  double pl = it->GetPersistenceLength();
  double dr = it->GetDriving();
  double nspec = GetNSpec();
  double theory;
  if (params_->n_dim == 2) {
    theory = l * pl * 4.0 - 8.0 * pl * pl * (1-exp(-0.5*l/pl));
  }
  else {
    theory = l * pl * 2.0 - 2.0 * pl * pl * (1-exp(-l/pl));
  }
  mse2e_file_ << l << " " << d << " " << cl << " " << pl << " " << dr << " " << params_->n_dim << " " << params_->n_steps << " " << nspec << " " << params_->delta << " " << theory << "\n";
  mse2e_file_ << "num_bead_springs_averaged mse2e_mean mse2e_std_err\n";
  mse2e_ = 0.0;
  mse2e2_ = 0.0;
  n_samples_ = 0;
}

void BeadSpringSpecies::InitThetaAnalysis() {
  // TODO Should check to make sure the same lengths, child lengths, persistence lengths, etc are used for each bead_spring in system.
  std::string fname = params_->run_name;
  fname.append("_bead_spring.theta");
  theta_file_.open(fname, std::ios::out);
  theta_file_ << "theta_analysis_file\n";
  theta_file_ << "length diameter bond_length persistence_length n_bead_springs n_bonds n_steps n_spec delta n_dim \n";
  double l, cl, pl, dr, d;
  int nbonds;
  int nmembers = members_.size();
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    l = it->GetLength();
    d = it->GetDiameter();
    cl = it->GetBondLength();
    pl = it->GetPersistenceLength();
    dr = it->GetDriving();
    nbonds = it->GetNBonds();
  }
  int nspec = GetNSpec();
  theta_file_ << l << " " << d << " " << cl << " " << pl << " " << dr << " " << nmembers << " " << nbonds << " " << params_->n_steps << " " << nspec << " " << params_->delta << " " << params_->n_dim << " " << "\n";
  theta_file_ << "cos_theta";
  for (int i=0; i<nbonds-1; ++i) {
    theta_file_ << " theta_" << i+1 << i+2;
  }
  theta_file_ << "\n";
  n_bins_ = 10000;
  int nfil = members_.size();
  theta_histogram_ = new int *[nbonds-1];
  for (int ibond=0; ibond<nbonds-1; ++ibond) {
    theta_histogram_[ibond] = new int[n_bins_];
    for (int ibin=0;ibin<n_bins_;++ibin) {
      theta_histogram_[ibond][ibin] = 0;
    }
  }
}

void BeadSpringSpecies::RunAnalysis() {
  // TODO Analyze conformation and ms end-to-end
  if (params_->bead_spring.theta_analysis) {
    if (params_->interaction_flag) {
      std::cout << "WARNING! Theta analysis running on interacting bead_springs!\n";
    }
    RunThetaAnalysis();
  }
  if (params_->bead_spring.lp_analysis) {
    RunMse2eAnalysis();
  }
  time_++;
}

void BeadSpringSpecies::RunMse2eAnalysis() {
  // Treat as though we have many spirals for now
  //if ( ! mse2e_file_.is_open()) {
    //early_exit = true;
    //std::cout << " Error! Problem opening file in RunMse2eAnalysis! Exiting.\n";
  //}
  //mse2e_file_ << time_;
  for (auto it=members_.begin(); it!= members_.end(); ++it) {
    double const * const head_pos = it->GetHeadPosition();
    double const * const tail_pos = it->GetTailPosition();
    double mse2e_temp = 0.0;
    for (int i=0; i<params_->n_dim; ++i) {
      double temp = (head_pos[i] - tail_pos[i]);
      mse2e_temp += temp*temp;
    }
    mse2e_ += mse2e_temp;
    mse2e2_ += mse2e_temp*mse2e_temp;
    //mse2e_file_ << " " << mse2e ;
  }
  //mse2e_ /= members_.size();
  //mse2e2_ /= members_.size();
  //mse2e_file_ << "\n";
  n_samples_++;
}


void BeadSpringSpecies::RunThetaAnalysis() {
  for (auto it=members_.begin(); it!=members_.end(); ++it) {
    std::vector<double> const * const thetas = it->GetThetas();
    for (int i=0; i<(it->GetNBonds()-1); ++i) {
      int bin_number = (int) floor( (1 + (*thetas)[i]) * (n_bins_/2) );
      if (bin_number == n_bins_) {
        bin_number = n_bins_-1;
      }
      else if (bin_number == -1) {
        bin_number = 0;
      }
      else if (bin_number > n_bins_ && bin_number < 0) {
        error_exit("Something went wrong in RunThetaAnalysis!");
      }
      theta_histogram_[i][bin_number]++;
    }
  }
}

void BeadSpringSpecies::FinalizeAnalysis() {
  if (spiral_file_.is_open()) {
    spiral_file_.close();
  }
  if (theta_file_.is_open()) {
    FinalizeThetaAnalysis();
    theta_file_.close();
  }
  if (mse2e_file_.is_open()) {
    FinalizeMse2eAnalysis();
    mse2e_file_.close();
  }
}

void BeadSpringSpecies::FinalizeMse2eAnalysis() {
  int num = members_.size();
  mse2e_file_ << num << " ";
  mse2e_ /= n_samples_*num;
  mse2e2_ /= n_samples_*num;
  mse2e_file_ << mse2e_ << " ";
  mse2e_file_ << sqrt((mse2e2_ - mse2e_*mse2e_)/(num*n_samples_)) << "\n";
}

void BeadSpringSpecies::FinalizeThetaAnalysis() {
  int nbonds = members_[members_.size()-1].GetNBonds();
  for (int i=0; i<n_bins_; ++i) {
    double axis = (2.0/n_bins_)*i - 1;
    theta_file_ << " " << axis;
    for (int ibond=0; ibond<nbonds-1; ++ibond) {
      theta_file_ << " " << theta_histogram_[ibond][i];
    }
    theta_file_ << "\n";
  }

  for (int ibond=0; ibond<nbonds-1; ++ibond) {
    delete[] theta_histogram_[ibond];
  }
  delete[] theta_histogram_;
}

