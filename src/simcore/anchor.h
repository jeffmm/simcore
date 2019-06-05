#ifndef _SIMCORE_ANCHOR_H_
#define _SIMCORE_ANCHOR_H_

#include "species.h"
#include "bond.h"
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

// Harmonic anchor structure
class Anchor : public Site {
  public:
    bool alignment_potential_,
         unbind_,
         bound_,
         walker_,
         diffuse_,
         active_;

    int step_direction_;
    double position_[3],
           orientation_[3],
           force_[3],
           torque_[3],
           k_spring_,
           k_align_,
           theta_, //euler angles
           phi_, //euler angles
           spring_length_,
           bond_length_,
           bond_lambda_,
           mesh_lambda_,
           velocity_,
           max_velocity_,
           diffusion_,
           f_spring_max_;
  public:
    Anchor();
    void Unbind() {
      unbind_ = true;
    }
    void Init();
    bool UpdatePriors();
    void UpdatePosition();
    void Diffuse();
    void DiffuseBound();
    void Activate();
    void Deactivate();
    void CheckNearBoundary();
    void CheckNearBuddingBoundary();
    void AnchorBoundary(double * anchor_point);
    void DetachBoundary();
    void ApplyAnchorForces();
    void AttachToBond(directed_bond, double lambda, double mesh_lambda);
    bool SwitchBonds(bool next_bond, double lambda);
    void UpdateAnchorPosition() {}
    void SetDiffusion();
    void SetWalker(int dir,double walk_v);
    void Walk();
    void AttachBondRandom(Bond * b, double mesh_lambda);
    double const GetMeshLambda();
};


typedef std::vector<Anchor>::iterator anchor_iterator;
typedef std::vector<std::pair<std::vector<Anchor>::iterator, std::vector<Anchor>::iterator> > anchor_chunk_vector;

class AnchorSpecies : public Species<Anchor> {
  protected:
    void CalculateBinding();
  public:
    AnchorSpecies() : Species() {
      SetSID(species_id::anchor);
    }
    void UpdatePositions() {
#ifdef ENABLE_OPENMP
      int max_threads = omp_get_max_threads();
      anchor_chunk_vector chunks;
      chunks.reserve(max_threads); 
      size_t chunk_size= members_.size() / max_threads;
      anchor_iterator cur_iter = members_.begin();
      for(int i = 0; i < max_threads - 1; ++i) {
        anchor_iterator last_iter = cur_iter;
        std::advance(cur_iter, chunk_size);
        chunks.push_back(std::make_pair(last_iter, cur_iter));
      }
      chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks)
      {
#pragma omp for 
        for(int i = 0; i < max_threads; ++i)
          for(auto it = chunks[i].first; it != chunks[i].second; ++it)
            it->UpdatePosition();
      }
#else
      for (anchor_iterator it=members_.begin(); it!=members_.end(); ++it) 
        it->UpdatePosition();
#endif
    }
};

#endif
