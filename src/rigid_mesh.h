#ifndef _SIMCORE_RIGID_MESH_H_
#define _SIMCORE_RIGID_MESH_H_

#include "mesh.h"

class RigidMesh : public Mesh {
  protected:
    double ix_position_[3],
           ix_prev_position_[3],
           ix_scaled_position_[3],
           ix_orientation_[3],
           ix_length_,
           ix_diameter_;
  public:
    RigidMesh();
    void SetInteractorLength(double const len);
    void SetInteractorDiameter(double const d);
    void SetInteractorPosition(double const * const pos);
    void SetInteractorPrevPosition(double const * const ppos);
    void SetInteractorScaledPosition(double const * const spos);
    void SetInteractorOrientation(double const * const u);
    virtual double const GetInteractorLength();
    virtual double const GetInteractorDiameter();
    virtual double const * const GetInteractorPosition();
    virtual double const * const GetInteractorPrevPosition();
    virtual double const * const GetInteractorScaledPosition();
    virtual double const * const GetInteractorOrientation();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
};

#endif
