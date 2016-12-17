#ifndef _SIMCORE_BR_BEAD_H_
#define _SIMCORE_BR_BEAD_H_

#include "object.h"
#include "auxiliary.h"
#include <yaml-cpp/yaml.h>
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

class BrBead : public Simple {
  private:
    double diffusion_,
           rand_sigma_rot_,
           gamma_rot_,
           driving_factor_,
           body_frame_[6];
  public:
    BrBead(system_parameters *params, space_struct *space, long seed, SID sid) : Simple(params, space, seed, sid) {
      diameter_=params->br_bead_diameter;
      driving_factor_ = params->driving_factor;
      SetDiffusion();
    }
    void SetDiffusion();
    void KickBead();
    void UpdatePosition();
    void GetBodyFrame();
    void Driving();
    void Translate();
    void Rotate();
    void Init();
};

#include "species.h"
class BrBeadSpecies : public Species<BrBead> {
  public:
    BrBeadSpecies() : Species() {
      SetSID(SID::br_bead);
    }
    BrBeadSpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::br_bead);
    }
    void Init() {Species::Init();}
    void UpdatePositions() {
#ifdef ENABLE_OPENMP
      int max_threads = omp_get_max_threads();
      std::vector<std::pair<std::vector<BrBead*>::iterator, std::vector<BrBead*>::iterator> > chunks;
      chunks.reserve(max_threads); 
      size_t chunk_size= members_.size() / max_threads;
      auto cur_iter = members_.begin();
      for(int i = 0; i < max_threads - 1; ++i) {
        auto last_iter = cur_iter;
        std::advance(cur_iter, chunk_size);
        chunks.push_back(std::make_pair(last_iter, cur_iter));
      }
      chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks)
      {
#pragma omp for 
        for(int i = 0; i < max_threads; ++i)
          for(auto it = chunks[i].first; it != chunks[i].second; ++it)
            (*it)->UpdatePosition();
      }
#else
      for (auto it=members_.begin(); it!=members_.end(); ++it) 
        (*it)->UpdatePosition();
#endif
    }

    void Configurator();
};

#endif // _SIMCORE_BR_BEAD_H_

