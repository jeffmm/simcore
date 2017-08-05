#ifndef _SIMCORE_OBJECT_H_
#define _SIMCORE_OBJECT_H_

#include "definitions.h"
#include "parameters.h"
#include "rng.h"

class Object {
  private:
    static int _next_oid_;
    static long _seed_;
  protected:
    static system_parameters const * params_;
    static int n_dim_;
    static double delta_;
    RNG rng_;
    int oid_;
    double position_[3],
           scaled_position_[3],
           prev_position_[3],
           orientation_[3],
           force_[3],
           torque_[3],
           length_,
           diameter_;
    graph_struct g_struct;
  public:
    Object();
    static void SetNDim(int n_dim) {n_dim_ = n_dim;}
    static void SetDelta(double delta) {delta_ = delta;}
    static void SetParameters(system_parameters * params) {params_=params;}
    static void SetSeed(long seed) {_seed_ = seed;}
    int GetOID();
    double const * const GetPosition();
    double const * const GetOrientation();
    double const GetLength();
    double const GetDiameter();
    void SetLength(double l);
    void SetDiameter(double d);
    void SetPosition(double * pos);
    void SetOrientation(double * u);
    void ZeroForce();
    virtual void UpdatePosition();
    virtual void UpdatePrevPosition();
    virtual void Report();
    // Main draw function, return struct of graphics info
    virtual void Draw(std::vector<graph_struct*> * graph_array);
};
#endif
