#ifndef _SIMCORE_SPHERE_LINE_WELL_H_
#define _SIMCORE_SPHERE_LINE_WELL_H_

#include "auxiliary.h"
#include "object.h"
#include "potential_base.h"

class SphereLineWell : public PotentialBase {
  protected:
    double rcut_touch_;
    double rcut_touch2_;
  public:
    SphereLineWell() : PotentialBase(nullptr, 0.0, 0.0) {
      pot_name_ = "SphereLineWell";
    }
    SphereLineWell(double pDepth, space_struct* pSpace, double pRcut, double pFcut) : PotentialBase(pSpace, pRcut, pFcut) {
      pot_name_ = "SphereLineWell";
    }
    virtual void Print() {
      PotentialBase::Print();
      std::cout << "\t{linepot: " << (4./3.)*M_PI*rcut_touch_*rcut_touch_*rcut_touch_ << "}\n";
    }

    virtual void CalcPotential(interactionmindist *idm,
                               Simple *part1,
                               Simple *part2,
                               double *fpote) {
      std::fill(fpote, fpote + n_dim_ + 1, 0.0);
      // 0 force, return depth when inside well
      double dr_mag2 = idm->dr_mag2;
      if (dr_mag2 < rcut_touch2_) {
        double a = sqrt(rcut_touch2_ - dr_mag2);
        double *rcontact;
        double rlength;
        // We have to know which one is the point, and which is the sphereo
        // cylinder
        if (part1->GetRigidLength() > 0.0) {
          rcontact = idm->contact1;
          rlength = part1->GetRigidLength();
        } else if (part2->GetRigidLength() > 0.0) {
          rcontact = idm->contact2;
          rlength = part2->GetRigidLength();
        } else {
          printf("ERROR, attempting point-line potential on non-point line\n");
          exit(1);
        }

        // Back calculate mu
        double rcontact2 = 0.0;
        for (int i = 0; i < n_dim_; ++i) {
          rcontact2 += SQR(rcontact[i]);
        }
        double mu = sqrt(rcontact2);
        double rmax = mu + a;
        if (rmax > 0.5 * rlength)
          rmax = 0.5 * rlength;
        else if (rmax < -0.5 * rlength)
          rmax = -0.5 * rlength;

        double rmin = mu - a;
        if (rmin > 0.5 * rlength)
          rmin = 0.5 * rlength;
        else if (rmin < -0.5 * rlength)
          rmin = -0.5 * rlength;

        fpote[n_dim_] = (rmax - rmin) / (4.0/3.0 * M_PI * rcut_touch_*rcut_touch_*rcut_touch_);
      } else {
        fpote[n_dim_] = 0.0;
      }
    }

    virtual void Init(space_struct *pSpace, int ipot, YAML::Node &node) {
      PotentialBase::Init(pSpace, ipot, node);

      // Now, let's look at the particular yaml node we are supposed to be interested in
      rcut_   = node["potentials"][ipot]["rcut"].as<double>();
      fcut_   = node["potentials"][ipot]["fcut"].as<double>();

      rcut2_ = rcut_*rcut_;
    }

    virtual void Init(space_struct *pSpace, YAML::Node *subnode) {
      YAML::Node node = *subnode;
      PotentialBase::Init(pSpace, &node);

      // Now, let's look at the particular yaml node we are supposed to be interested in
      rcut_       = node["rcut"].as<double>();
      fcut_       = node["fcut"].as<double>();
      rcut_touch_ = node["rcut_touch"].as<double>();

      rcut_touch2_ = rcut_touch_*rcut_touch_;
      rcut2_ = rcut_*rcut_;
    }

};

#endif

