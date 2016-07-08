#ifndef _SIMCORE_INTERACTION_ENGINE_H_
#define _SIMCORE_INTERACTION_ENGINE_H_

#include <unordered_map>
#include <memory>

#include "auxiliary.h"
#include "species.h"
#include "minimum_distance.h"
#include "neighbor_list_generic.h"
#include "particle_tracking.h"
#include "potential_manager.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class InteractionEngine {
  public:
    InteractionEngine() {}
    ~InteractionEngine() {
      delete[] frc_;
      delete[] trqc_;
      delete[] prc_energy_;
      delete[] kmc_energy_;
      oid_position_map_.clear();
    }

    void Init(space_struct *pSpace, ParticleTracking *pTracking, double pSkin);
    void InitPotentials(PotentialManager *pPotentials);
    void InitMP();

    void Interact();
    void InteractParticlesMP(int &idx,
                             int &jdx,
                             double **fr,
                             double **tr,
                             double *pr_energy,
                             double *kmc_energy);
    void KMCParticlesMP(neighbor_t *neighbor, int &idx, int &jdx);
    /*void InteractParticlesMP(neighbor_t *neighbor,
                             Simple *part1,
                             Simple *part2,
                             double **fr,
                             double **tr,
                             double *pr_energy,
                             double *kmc_energy);*/
    void ReduceParticlesMP();
    void Print();
    void Dump();

  private:
    int ndim_;
    int nperiodic_;
    int nthreads_;
    int nsimples_;

    double max_rcut_;
    double skin_;
    double box_[3];

    double* frc_ = nullptr; // Force superarray
    double* trqc_ = nullptr; // torque superarray
    double* prc_energy_ = nullptr; // Potential energy superarray
    double* kmc_energy_ = nullptr; // kmc energy superarray
    std::unordered_map<int, int> oid_position_map_; // oid to position mapping!!!

    space_struct *space_;
    ParticleTracking *tracking_;
    PotentialManager *potentials_;
    std::vector<Simple*>* simples_;
};

// Helper functions for the setup/teardown of the MP
// regions
namespace fmmph {
  // Init the MP region (all pointers)
  // HAve to main ref to pointers(gross)
  inline void InitMPRegion(int &tid,
                           int &nparticles,
                           double **frc,
                           double **trqc,
                           double **prc_energy,
                           double **kmc_energy,
                           double ***fr,
                           double ***tr,
                           double **pr,
                           double **kmcr) {
    // Set up the pointers to the force and torque superarrays
    // Do all 3 dimensions to be safe
    // Do 1d first, it's easier
    (*pr)   = (*prc_energy) + (tid * nparticles);
    (*kmcr) = (*kmc_energy) + (tid * nparticles);
    for (int i = 0; i < nparticles; ++i) {
      (*pr)[i]    = 0.0;
      (*kmcr)[i]  = 0.0;
    }

    // Multidimensional now
    for (int idim = 0; idim < 3; ++idim) {
      (*fr)[idim] = (*frc) + ((3*tid+idim)*nparticles);
      (*tr)[idim] = (*trqc)+ ((3*tid+idim)*nparticles);
    }
    for (int idx = 0; idx < nparticles; ++idx) {
      for (int idim = 0; idim < 3; ++idim) {
        (*fr)[idim][idx] = 0.0;
        (*tr)[idim][idx] = 0.0;
      }
    }
  }

  // Reduce the openmp region to the superarrays
  inline void ReduceMPRegion(int &tid,
                             int &nparticles,
                             int &nthreads,
                             double **frc,
                             double **trqc,
                             double **prc_energy,
                             double **kmc_energy) {
    int i = 1 + (3 * nparticles / nthreads);
    int ii = 1 + (nparticles / nthreads);

    int fromidx = tid * i;
    int fromidx2 = tid * ii;

    int toidx = fromidx + i;
    int toidx2 = fromidx2 + ii;

    if (toidx > 3*nparticles) toidx = 3*nparticles;
    if (toidx2 > nparticles) toidx2 = nparticles;

    // Reduce the forces
    for (i = 1; i < nthreads; ++i) {
      int offs = 3 * i * nparticles;

      for (int j = fromidx; j < toidx; ++j) {
        (*frc)[j] += (*frc)[offs+j];
        (*trqc)[j] += (*trqc)[offs+j];
      }
    }

    // Reduce the energies
    for (ii = 1; ii < nthreads; ++ii) {
      int offs = ii * nparticles;

      for (int jj = fromidx2; jj < toidx2; ++jj) {
        (*prc_energy)[jj] += (*prc_energy)[offs+jj];
        (*kmc_energy)[jj] += (*kmc_energy)[offs+jj];
      }
    }
  }
}

#endif
