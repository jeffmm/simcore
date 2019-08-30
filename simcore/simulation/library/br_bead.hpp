#ifndef _SIMCORE_BR_BEAD_H_
#define _SIMCORE_BR_BEAD_H_

#include "species.hpp"
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

class BrBead : public Object {
 protected:
  bool stoch_flag_;
  double gamma_trans_, gamma_rot_, diffusion_, driving_factor_;
  void ApplyForcesTorques();
  void ApplyBoundaryForces();
  void InsertBrBead();
  void SetDiffusion();
  void Translate();
  void Rotate();
  void Integrate();

 public:
  BrBead();
  void Init();
  void UpdatePosition();
  virtual void GetInteractors(std::vector<Object *> *ix);
  virtual int GetCount();
  virtual void Draw(std::vector<graph_struct *> *graph_array);
  virtual void ZeroForce();
};

typedef std::vector<BrBead>::iterator br_bead_iterator;
typedef std::vector<
    std::pair<std::vector<BrBead>::iterator, std::vector<BrBead>::iterator> >
    br_bead_chunk_vector;

class BrBeadSpecies : public Species<BrBead> {
 public:
  BrBeadSpecies() : Species() { SetSID(species_id::br_bead); }
  void Init(system_parameters *params, space_struct *space, long seed) {
    Species::Init(params, space, seed);
    sparams_ = &(params_->br_bead);
  }
  void UpdatePositions() {
#ifdef ENABLE_OPENMP
    int max_threads = omp_get_max_threads();
    br_bead_chunk_vector chunks;
    chunks.reserve(max_threads);
    size_t chunk_size = members_.size() / max_threads;
    br_bead_iterator cur_iter = members_.begin();
    for (int i = 0; i < max_threads - 1; ++i) {
      br_bead_iterator last_iter = cur_iter;
      std::advance(cur_iter, chunk_size);
      chunks.push_back(std::make_pair(last_iter, cur_iter));
    }
    chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks)
    {
#pragma omp for
      for (int i = 0; i < max_threads; ++i)
        for (auto it = chunks[i].first; it != chunks[i].second; ++it)
          it->UpdatePosition();
    }
#else
    for (br_bead_iterator it = members_.begin(); it != members_.end(); ++it)
      it->UpdatePosition();
#endif
  }
};

#endif  // _SIMCORE_BR_BEAD_H_
