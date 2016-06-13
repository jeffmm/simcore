#ifndef _SIMCORE_SPACE_PROPERTIES_H_
#define _SIMCORE_SPACE_PROPERTIES_H_

#include <gsl/gsl_rng.h>
#include <math.h>
#include "auxiliary.h"
#include "macros.h"

typedef enum {
  SPHERE = 0,
  BOX = 1,
  SNOWMAN = 2,
} boundary_type_t;

class SpaceProperties {
  private:
    boundary_type_t boundary_type_;
    int n_dim_,
        n_periodic_;
    double radius_,
           d_radius_, // radius of daughter cell, if applicable
           m_d_dist_, // mother-daughter cell separation
           intersect_height_, // height of plane at mother-daughter cell intersection
           intersect_radius_, // radius of circle defined at mother-daughter cell intersection
           r_cutoff_, // WCA potential cutoff at boundary
           volume_,
           uc_volume_,
           v_ratio_,
           **a_, // direct lattice vector
           **b_, // reciprocal lattice vector
           *a_perp_, // perp dist between opposite unit cell faces
           **uc_,
           **uc_inv_; // inverse unit cell matrix
    rng_properties rng_;
    system_parameters *params_;
    space_struct s_struct;

  public:
    SpaceProperties();
    void Init(system_parameters *params, long seed);
    void InitUnitCell();
    void ClearUnitCell();
    void Clear();
    void CalculateVolume();
    void RandomCoordinate(double *vec);
    void RandomCoordinate(double *vec, double buffer);
    int const GetDim();
    int const GetPeriodic();
    double const GetIntersectHeight();
    double const GetIntersectRadius();
    double const GetRadius();
    double const GetDRadius();
    double const GetMDDist();
    double const GetVolume();
    double const * const * const GetA();
    double const * const * const GetB();
    double const * const GetAPerp();
    double const * const * const GetUnitCell();
    double const * const * const GetUnitCellInv();
    boundary_type_t GetType();
    std::string GetTypeString();
    bool CheckInBounds(double *vec, double buffer);
    bool CheckSegmentInBounds(double *vec1, double *vec2, double buffer);
    space_struct * GetStruct();
    void UpdateSpaceStruct();
    //static void InsertRandom(double pos[3], double buffer);
};

#endif // _SIMCORE_SPACE_PROPERTIES_H_
