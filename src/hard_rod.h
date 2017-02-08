#ifndef _SIMCORE_HARD_ROD_H_
#define _SIMCORE_HARD_ROD_H_

#include "site.h"
#include "bond.h"
#include "species.h"
#include "auxiliary.h"
#include <yaml-cpp/yaml.h>
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

class HardRod : public Composite<Site,Bond> {

  private:
    int n_bonds_;
    double max_length_,
           min_length_,
           max_child_length_,
           child_length_,
           friction_par_,
           friction_perp_,
           friction_rot_,
           rand_sigma_par_,
           rand_sigma_perp_,
           rand_sigma_rot_,
           body_frame_[6],
           driving_factor_;
    poly_state_t poly_state_;
    void UpdateSitePositions();
    void UpdateBondPositions();
    void SetDiffusion();
    void UpdateOrientation();
    void GetBodyFrame();
    void AddRandomDisplacement();
    void ApplyForcesTorques();

  public:
    HardRod(system_parameters *params, space_struct * space, long seed, SID sid) 
      : Composite(params, space, seed, sid) {
        length_ = params->hard_rod.length;
        diameter_ = params->hard_rod.diameter;
        max_length_ = params->hard_rod.max_length;
        min_length_ = params->hard_rod.min_length;
        max_child_length_ = params->hard_rod.max_child_length;
        driving_factor_ = params->hard_rod.driving_factor;
        // Initialize end sites
        for (int i=0; i<2; ++i) {
          Site s(params, space, gsl_rng_get(rng_.r), GetSID());
          s.SetCID(GetCID());
          elements_.push_back(s);
        }
        // Initialize bonds
        n_bonds_ = (int) ceil(length_/max_child_length_);
        child_length_ = length_/n_bonds_;
        for (int i=0; i<n_bonds_; ++i) {
          Bond b(params, space, gsl_rng_get(rng_.r), GetSID());
          b.SetCID(GetCID());
          b.SetRID(GetRID());
          v_elements_.push_back(b);
        }
      }
    virtual void Init();
    virtual void Integrate();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePosition();
    virtual void Dump();
    void ScalePosition();

};

class HardRodSpecies : public Species<HardRod> {
  protected:
    //void InitPotentials(system_parameters *params);
    double max_length_;
    double min_length_;
  public:
    HardRodSpecies() : Species() {
      SetSID(SID::hard_rod);
    }
    HardRodSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      SetSID(SID::hard_rod);
      //InitPotentials(params);
      max_length_ = params_->hard_rod.max_length;
      min_length_ = params_->hard_rod.min_length;
    }
    double const GetMaxLength() {return max_length_;}
    double const GetMinLength() {return min_length_;}

    void UpdatePositions() {
#ifdef ENABLE_OPENMP
      int max_threads = omp_get_max_threads();
      std::vector<std::pair<std::vector<HardRod*>::iterator, std::vector<HardRod*>::iterator> > chunks;
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
      for (auto it=members_.begin(); it != members_.end(); ++it) 
        (*it)->UpdatePosition();
#endif
    }
    // Special insertion routine
    void Configurator();
};

#endif // _SIMCORE_HARD_ROD_H_
