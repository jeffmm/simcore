#ifndef _SIMCORE_INTERACTION_ENGINE_H_
#define _SIMCORE_INTERACTION_ENGINE_H_

#include <unordered_map>
#include <memory>

#include "anchor_list_generic.h"
#include "auxiliary.h"
#include "interaction.h"
#include "particle_engine.h"
#include "species.h"
#include "minimum_distance.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class InteractionEngine {
  public:
    InteractionEngine() {}
    virtual ~InteractionEngine() {}

    void Dump();
    void Init(space_struct *pSpace,
              ParticleEngine *pTrackEngine,
              std::vector<interaction_t> *pInteractions);
    void InitMP();
    void Interact();

  protected:
    int ndim_;
    int nperiodic_;
    int nthreads_;
    int nsimples_;
    int nspecies_;

    double box_[3];

    // Superarrays
    double *frc_  = nullptr;
    double *trqc_ = nullptr;
    double *prc_energy_ = nullptr;
    double *virial_ = nullptr;

    space_struct *space_;
    ParticleEngine *ptrack_;
    std::vector<interaction_t> *interactions_;
    std::vector<Simple*> *simples_;
    std::vector<SpeciesBase*> *species_;
    std::unordered_map<int, int> *oid_position_map_;
    al_set *anchors_;

    std::map<SID, int> spec_ind_map_;

    void AttachParticleEngine();
    void DumpInteractions();
    void DumpSpecies();
    void ReduceParticlesMP();

    // All interaction types
    void InteractParticlesExternalMP(interaction_t **pix,
                                     double **fr,
                                     double **tr,
                                     double *pe,
                                     double **virial);
    void InteractParticlesKMCMP(interaction_t **pix,
                                double **fr,
                                double **tr,
                                double *pe,
                                double **virial);
    void InteractParticlesInternalMP(interaction_t **pix,
                                double **fr,
                                double **tr,
                                double *pe);
    void InteractParticlesBoundaryMP(interaction_t **pix,
                                     double **fr,
                                     double **tr,
                                     double *pe,
                                     double **virial);
    void InteractParticlesTetherMP(interaction_t **pix,
                                   double **fr,
                                   double **tr,
                                   double *pe,
                                   double **virial);



};

// Helper functions for the setup/teardown of the MP
// regions
namespace ieh {
  // Init the MP region (all pointers)
  // HAve to main ref to pointers(gross)
  inline void InitMPRegion(int &tid,
                           int &nparticles,
                           int &nspecies,
                           double **frc,
                           double **trqc,
                           double **prc_energy,
                           double **virialc,
                           double ***fr,
                           double ***tr,
                           double **pr,
                           double ***virial) {
    // Set up the pointers to the force and torque superarrays
    // Do all 3 dimensions to be safe
    // Do 1d first, it's easier
    (*pr)   = (*prc_energy) + (tid * nparticles);
    for (int i = 0; i < nparticles; ++i) {
      (*pr)[i]    = 0.0;
    }

    // Multidimensional now
    for (int idim = 0; idim < 3; ++idim) {
      (*fr)[idim] = (*frc) + ((3*tid+idim)*nparticles);
      (*tr)[idim] = (*trqc)+ ((3*tid+idim)*nparticles);
      for (int jdim = 0; jdim < 3; ++jdim) 
        (*virial)[3*idim+jdim] = (*virialc)+((9*tid+3*idim+jdim)*nspecies);
    }
    for (int idx = 0; idx < nparticles; ++idx) {
      for (int idim = 0; idim < 3; ++idim) {
        (*fr)[idim][idx] = 0.0;
        (*tr)[idim][idx] = 0.0;
      }
    }
    for (int ids = 0; ids < nspecies; ++ids) {
      for (int idim = 0; idim < 3; ++idim) 
        for (int jdim = idim; jdim < 3; ++jdim) 
          (*virial)[3*idim+jdim][ids] = (*virial)[3*jdim+idim][ids] = 0.0;
    }
  }

  // Reduce the openmp region to the superarrays
  inline void ReduceMPRegion(int &tid,
                             int &nparticles,
                             int &nspecies,
                             int &nthreads,
                             double **frc,
                             double **trqc,
                             double **virialc,
                             double **prc_energy) {
    int i = 1 + (3 * nparticles / nthreads);
    int ii = 1 + (nparticles / nthreads);
    int is = 1 + (9 * nspecies / nthreads);

    int fromidx = tid * i;
    int fromidx2 = tid * ii;
    int fromids = tid * is;

    int toidx = fromidx + i;
    int toidx2 = fromidx2 + ii;
    int toids = fromids + is;

    if (toidx > 3*nparticles) toidx = 3*nparticles;
    if (toidx2 > nparticles) toidx2 = nparticles;
    if (toids > 9*nspecies) toids = 9*nspecies;

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
      }
    }

    for (is = 1; is < nthreads; ++is){
      int offs = is * nspecies;
      for (int js = fromids; js < toids; ++js) 
        (*virialc)[js] += (*virialc)[offs+js];
    }
  }
}

#endif /* _SIMCORE_INTERACTION_ENGINE_H_ */
