#ifndef _SIMCORE_OBJECT_H_
#define _SIMCORE_OBJECT_H_

#include "auxiliary.h"

class Object {

  private:
    unsigned int oid_;
    static unsigned int next_oid_;
        
  protected:
    int n_dim_;
    double *position_,
           *orientation_,
           *velocity_,
           *force_,
           radius_,
           length_;
    std::vector<interaction> interactions;

  public:
    Object();
    Object(int n_dim);
    virtual ~Object();
    const unsigned int GetObjectID() const;
    virtual void ZeroForce();
    virtual void SetPosition(double *new_position);
    virtual void SetOrientation(double *new_orientation);
    virtual void SetVelocity(double *new_velocity);
    virtual void ApplyInteractions();
    void Init(int n_dim);
    void AddInteraction(interaction new_interaction);
    double *GetPosition();
    double *GetOrientation();
    double *GetVelocity();
    double *GetForce();
    virtual void SetRadius(double new_radius);
    virtual void SetLength(double new_length);
    double GetRadius() const;
    double GetLength() const;
};

#endif // _SIMCORE_OBJECT_H_
