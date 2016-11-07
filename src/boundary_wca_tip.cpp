#include "boundary_wca_tip.h"

#include "object.h"

void BoundaryWCATip::Print() {
  PotentialBase::Print();
  std::cout << std::setprecision(16)
    << "\t{eps: " << eps_ << "}\n"
    << "\t{sigma: " << sigma_ << "}\n"
    << "\t{c6: " << c6_ << "}\n"
    << "\t{c12: " << c12_ << "}\n"
    << "\t{conf_radius: " << conf_radius_ << "}\n"
    << "\t{tip: " << tip_ << "}"
    << std::endl;
}

void BoundaryWCATip::CalcPotential(interactionmindist *idm,
                                Simple *part1,
                                Simple *part2,
                                double *fpote) {
  std::fill(fpote, fpote + n_dim_ + 1, 0.0);
  // Do slightly differently, since we have only one particle,
  // and no distance calculation
  // Calculate dr to the boundary
  // Also, assign rcontact1 to idm for use later...
  double rmag2 = 0.0;
  double rx[3] = {0.0, 0.0, 0.0};
  double ux[3] = {0.0, 0.0, 0.0};
  double l = 0.0;
  double rcontact[3] = {0.0, 0.0, 0.0};
  double prefactor = 0.0;
  if (tip_ == 0)
    prefactor = -0.5;
  else
    prefactor = 0.5;
  std::copy(part1->GetRigidPosition(), part1->GetRigidPosition()+3, rx);
  std::copy(part1->GetRigidOrientation(), part1->GetRigidOrientation()+3, ux);
  l = part1->GetRigidLength();
  for (int i = 0; i < n_dim_; ++i) {
    rcontact[i] = prefactor * l * ux[i];
    idm->contact1[i] = rcontact[i];
    rmag2 += SQR(rx[i] + rcontact[i]);
  }

  if (rmag2 > conf_radius2_) {
    // Outside the boundary, interact!
    double rmag = sqrt(rmag2);

    double delta_r = -rmag + conf_radius_ + rcut_;
    double rinv = 1.0/(delta_r*delta_r);
    double r6 = rinv*rinv*rinv;

    double ffac = (12.0*c12_*r6 - 6.0*c6_)*r6*rinv;
    // Cut off the force at fcut
    if (ABS(ffac) > fcut_) {
      ffac = SIGNOF(ffac) * fcut_;
    }

    for (int i = 0; i < n_dim_; ++i) {
      fpote[i] = -ffac * delta_r * (rx[i] + rcontact[i]) / rmag;
    }

    fpote[n_dim_] = r6*(c12_*r6 - c6_) + eps_;
  }
}

void BoundaryWCATip::Init(space_struct *pSpace, int ipot, YAML::Node &node) {
    PotentialBase::Init(pSpace, ipot, node);

    // Now, let's look at the particular yaml node we are supposed to be interested in
    eps_    = node["potentials"][ipot]["eps"].as<double>();
    sigma_  = node["potentials"][ipot]["sigma"].as<double>();
    fcut_   = node["potentials"][ipot]["fcut"].as<double>();
    tip_    = node["potentials"][ipot]["tip"].as<int>();

    // For WCA potentials, the rcutoff is actually important, as it must be
    // restricted to be at 2^(1/6)sigma

    rcut_ = pow(2.0, 1.0/6.0)*sigma_;

    rcut2_ = rcut_*rcut_;
    c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
    c6_  = 4.0 * eps_ * pow(sigma_,  6.0);

    // The radius is set by the system boundary!!!!
    // Taken from the definition in bob
    conf_radius_ = 0.5 * space_->unit_cell[0][0] - rcut_ + 0.5;
    conf_radius2_ = conf_radius_ * conf_radius_;
}

void BoundaryWCATip::Init(space_struct *pSpace, YAML::Node *subnode) {
  YAML::Node node = *subnode;
  PotentialBase::Init(pSpace, &node);

  // Now, let's look at the particular yaml node we are supposed to be interested in
  eps_    = node["eps"].as<double>();
  sigma_  = node["sigma"].as<double>();
  fcut_   = node["fcut"].as<double>();
  tip_    = node["tip"].as<int>();

  // For WCA potentials, the rcutoff is actually important, as it must be
  // restricted to be at 2^(1/6)sigma

  rcut_ = pow(2.0, 1.0/6.0)*sigma_;

  rcut2_ = rcut_*rcut_;
  c12_ = 4.0 * eps_ * pow(sigma_, 12.0);
  c6_  = 4.0 * eps_ * pow(sigma_,  6.0);

  // The radius is set by the system boundary!!!!
  // Taken from the definition in bob
  conf_radius_ = 0.5 * space_->unit_cell[0][0] - rcut_ + 0.5;
  conf_radius2_ = conf_radius_ * conf_radius_;
}
