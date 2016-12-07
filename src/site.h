#ifndef _SIMCORE_SITE_H_
#define _SIMCORE_SITE_H_

#include "species.h"
#include "object.h"
#include "auxiliary.h"

class Site : public Simple {
  private:
    double tangent_[3];
    double random_force_[3];
  public:
    Site(system_parameters *params, space_struct *space, 
        long seed, SID sid) : Simple(params, space, seed, sid) {}
    void Init();
    void Draw(std::vector<graph_struct*> * graph_array);
    double const * const GetTangent() {return tangent_;}
    void SetTangent(double const * const tan) {std::copy(tan, tan+n_dim_, tangent_);}
    void SetRandomForce(double const * const f_rand) {std::copy(f_rand, f_rand+n_dim_, random_force_);}
    void AddRandomForce() {for (int i=0; i<n_dim_; ++i) force_[i] += random_force_[i];}
    double const * const GetRandomForce() {return random_force_;}
};

#endif // _SIMCORE_SITE_H_
