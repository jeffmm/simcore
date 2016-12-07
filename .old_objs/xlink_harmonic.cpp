#include "xlink_harmonic.h"

// Include xlinks...
#include "xlink_head.h"

void XlinkHarmonic::Print() {
  PotentialBase::Print();
  std::cout << " \tk:                  " << std::setprecision(16) << k_ << std::endl;
  std::cout << " \tequilibrium_length: " << std::setprecision(16) << r_equil_ << std::endl;
}

void XlinkHarmonic::CalcPotential(interactionmindist *idm,
                           Simple *part1,
                           Simple *part2,
                           double *fpote) {
  std::fill(fpote, fpote + n_dim_ + 1, 0.0);
  // Both heads are xlink_heads, cast them so we can check bound state
  XlinkHead *head0 = dynamic_cast<XlinkHead*>(part1);
  XlinkHead *head1 = dynamic_cast<XlinkHead*>(part2);
  auto boundstate = head0->GetBound();
  if (!head0->GetBound() || !head1->GetBound()) {
    return;
  }
  double rmag = sqrt(idm->dr_mag2);
  double *dr = idm->dr;

  double k = k_;
  double ffac;
  if (r_equil_ == 0.0) {
    ffac = k;
  } else {
    ffac = k * (1.0 - r_equil_ / sqrt(dot_product(n_dim_, dr, dr)));
  }
  double u = 0.0;

  if (ABS(ffac) > fcut_)
    ffac = SIGNOF(ffac) * fcut_;
  for (int i = 0; i < n_dim_; ++i)  {
    fpote[i] = dr[i]*ffac;
    u += 0.5 * SQR(fpote[i]);
  }
  u *= 0.5 / k;
  fpote[n_dim_] = u;
}

void XlinkHarmonic::Init(space_struct *pSpace, int ipot, YAML::Node &node) {
  PotentialBase::Init(pSpace, ipot, node);

  // Now, let's look at the particular yaml node we are supposed to be interested in
  rcut_     = node["potentials"][ipot]["rcut"].as<double>();
  k_        = node["potentials"][ipot]["k"].as<double>();
  r_equil_  = node["potentials"][ipot]["equilibrium_length"].as<double>();
  fcut_     = node["potentials"][ipot]["fcut"].as<double>();

  rcut2_ = rcut_*rcut_;
}

void XlinkHarmonic::Init(space_struct *pSpace, YAML::Node *subnode) {
  YAML::Node node = *subnode;
  PotentialBase::Init(pSpace, &node);

  // Now, let's look at the particular yaml node we are supposed to be interested in
  rcut_     = node["rcut"].as<double>();
  k_        = node["k"].as<double>();
  r_equil_  = node["equilibrium_length"].as<double>();
  fcut_     = node["fcut"].as<double>();

  rcut2_ = rcut_*rcut_;
}


