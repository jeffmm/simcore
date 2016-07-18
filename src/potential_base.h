#ifndef _SIMCORE_POTENTIAL_BASE_H_
#define _SIMCORE_POTENTIAL_BASE_H_

#include "auxiliary.h"
#include <yaml-cpp/yaml.h>

// Forward declare Simples
class Simple;

// Potential base class for all external(or maybe even internal) potentials
class PotentialBase {
  protected:
    bool can_overlap_ = false;
    bool is_kmc_ = false; // Is this a kinetic monte carlo potential (no forces?)
    int n_dim_;
    double rcut_, rcut2_; // Cutoff radius
    double fcut_; // Force cutoff
    SID kmc_target_ = SID::none; // KMC target (if there is one)
    std::string pot_name_;
    space_struct *space_;
  public:
    PotentialBase(space_struct* pSpace, double pRcut, double pFcut) : rcut_(pRcut), space_(pSpace), fcut_(pFcut) {
      pot_name_ = "Untitled";
      rcut2_ = rcut_ * rcut_;
      if (space_ != nullptr)
        n_dim_ = space_->n_dim;
      else
        n_dim_ = 0;
    }
    virtual ~PotentialBase() {}
    virtual void CalcPotential(interactionmindist *idm,
                               Simple *part1,
                               Simple *part2,
                               double *fpote) {
      printf("ERROR, CalcPotential not set correctly\n");
      exit(1);
    }
    double GetRCut() {return rcut_;}
    double GetRCut2() {return rcut2_;}
    const bool CanOverlap() { return can_overlap_; }
    const bool IsKMC() { return is_kmc_; }
    const SID GetKMCTarget() { return kmc_target_; }
    virtual void Print() {
      std::cout << pot_name_ << "\n";
      printf("\t{rcut:%2.2f}, {kmc: %s}, {kmc_target: %d}\n", rcut_, is_kmc_ ? "true" : "false",
            (int)kmc_target_);
    }

    virtual void Init(space_struct *pSpace, int ipot, YAML::Node &node) {
      space_ = pSpace;
      n_dim_ = space_->n_dim;

      if (node["potentials"][ipot]["overlap"]) {
        can_overlap_ = node["potentials"][ipot]["overlap"].as<bool>();
      }
      if (node["potentials"][ipot]["kmc"]) {
        is_kmc_ = node["potentials"][ipot]["kmc"].as<bool>();
        // Must specifcy a kmc_target for particles
        std::string kmc_sid_str = node["potentials"][ipot]["kmc_target"].as<std::string>();
        SID kmc_sid = StringToSID(kmc_sid_str);
        kmc_target_ = kmc_sid;
      }
  }
};

typedef std::pair<sid_pair, PotentialBase*> potential_pair;
typedef std::map<sid_pair, PotentialBase*> potential_map;

#endif // _SIMCORE_POTENTIAL_BASE_H_
