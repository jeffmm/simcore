#include "harmonic_tip.h"

#include "object.h"

void HarmonicTip::Print() {
  PotentialBase::Print();
  std::cout << " \tk:                  " << std::setprecision(16) << k_ << std::endl;
  std::cout << " \tequilibrium_length: " << std::setprecision(16) << r_equil_ << std::endl;
  std::cout << " \ttip:                " << std::setprecision(16) << tip_ << std::endl;
}

void HarmonicTip::CalcPotential(interactionmindist *idm,
                           Simple *part1,
                           Simple *part2,
                           double *fpote) {
  std::fill(fpote, fpote + n_dim_ + 1, 0.0);
  // Do slightly differently, since we have only one particle,
  // and no distance calculation
  // Calculate dr to the boundary
  // Also, assign rcontact1 to idm for use later...
  // part1 is the spindle pole body, part2 is the rod (tip)
  double rmag2 = 0.0;
  double rx[3] = {0.0, 0.0, 0.0};
  double sx[3] = {0.0, 0.0, 0.0};
  double ux[3] = {0.0, 0.0, 0.0};
  double l = 0.0;
  double rcontact[3] = {0.0, 0.0, 0.0};
  double prefactor = 0.0;
  if (tip_ == 0)
    prefactor = -0.5;
  else
    prefactor = 0.5;
  std::copy(part2->GetRigidPosition(), part2->GetRigidPosition()+3, rx);
  std::copy(part2->GetScaledPosition(), part2->GetScaledPosition()+3, rx);
  std::copy(part2->GetRigidOrientation(), part2->GetRigidOrientation()+3, ux);
  l = part1->GetRigidLength();
  // Get the position of the tip
  for (int i = 0; i < n_dim_; ++i) {
    rcontact[i] = prefactor * l * ux[i];
    rx[i] = rx[i] + rcontact[i];
  }
  double dr[3] = {0.0, 0.0, 0.0};
  separation_vector(n_dim_, space_->n_periodic, part1->GetRigidPosition(), part1->GetScaledPosition(),
      rx, sx, space_->unit_cell, dr);

  

  std::cout << "Harmonic Tip DEBUG, exiting\n";
  exit(1);




  // Slightly different, since we're using the tip of a rigid object
  //std::fill(fpote, fpote + n_dim_ + 1, 0.0);



  //double rmag = idm->dr_mag;
  //double *dr = idm->dr;

  //double k = k_;
  //double ffac;
  //if (r_equil_ == 0.0) {
  //  ffac = k;
  //} else {
  //  ffac = k * (1.0 - r_equil_ / sqrt(dot_product(n_dim_, dr, dr)));
  //}
  //double u = 0.0;

  //if (ABS(ffac) > fcut_)
  //  ffac = SIGNOF(ffac) * fcut_;
  //for (int i = 0; i < n_dim_; ++i)  {
  //  fpote[i] = dr[i]*ffac;
  //  u += 0.5 * SQR(fpote[i]);
  //}
  //u *= 0.5 / k;
  //fpote[n_dim_] = u;
}

void HarmonicTip::Init(space_struct *pSpace, int ipot, YAML::Node &node) {
  PotentialBase::Init(pSpace, ipot, node);

  // Now, let's look at the particular yaml node we are supposed to be interested in
  rcut_     = node["potentials"][ipot]["rcut"].as<double>();
  k_        = node["potentials"][ipot]["k"].as<double>();
  r_equil_  = node["potentials"][ipot]["equilibrium_length"].as<double>();
  fcut_     = node["potentials"][ipot]["fcut"].as<double>();
  tip_      = node["potentials"][ipot]["tip"].as<int>();

  rcut2_ = rcut_*rcut_;
}

