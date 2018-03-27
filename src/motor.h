#include "species.h"
#include "bond.h"
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

class Motor : public Site {

  protected:
    bool bound_,
         walker_,
         diffuse_,
         active_,
         anchored_;
    int step_direction_;
    double k_on_,
           k_off_,
           bond_length_,
           bond_lambda_,
           mesh_lambda_,
           velocity_,
           max_velocity_,
           diffusion_,
           f_spring_max_;
    Anchor anchor_;
  public:
    Motor();
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
    void UpdateMotorPosition() {}
    void SetDiffusion();
    void SetWalker(int dir,double walk_v);
    void Walk();
    void AttachBondRandom(Bond * b, double mesh_lambda);
    double const GetMeshLambda();
};

typedef std::vector<Motor>::iterator motor_iterator;
typedef std::vector<std::pair<std::vector<Motor>::iterator, std::vector<Motor>::iterator> > motor_chunk_vector;

class MotorSpecies : public Species<Motor> {
  protected:
    void CalculateBinding();
  public:
    void UpdatePositions() {
#ifdef ENABLE_OPENMP
      int max_threads = omp_get_max_threads();
      motor_chunk_vector chunks;
      chunks.reserve(max_threads); 
      size_t chunk_size= members_.size() / max_threads;
      motor_iterator cur_iter = members_.begin();
      for(int i = 0; i < max_threads - 1; ++i) {
        motor_iterator last_iter = cur_iter;
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
      for (motor_iterator it=members_.begin(); it!=members_.end(); ++it) 
        it->UpdatePosition();
#endif
    }
};

