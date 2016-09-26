#ifndef _SIMCORE_SPINDLE_POLE_BODY_H_
#define _SIMCORE_SPINDLE_POLE_BODY_H_

#include "anchor_list_generic.h"
#include "auxiliary.h"
#include "object.h"

class SpindlePoleBody : public Simple {

  protected:
    // Local frame of reference for SPB
    double u_anchor_[3]; // u toward center of nucleus, (-ranchor/R)
    double v_anchor_[3];
    double w_anchor_[3];
    double conf_rad_;

    double tau_local_[3]; // local torque about SPB axis
    double gamma_tra_;
    double gamma_rot_;
    double attach_diameter_;

  public:
    SpindlePoleBody(system_parameters *params, space_struct *space, long seed, SID sid) : Simple(params, space, seed, sid) {}
    ~SpindlePoleBody() {}
    SpindlePoleBody(const SpindlePoleBody& that) : Simple(that) {}
    SpindlePoleBody& operator=(SpindlePoleBody const& that) {Simple::operator=(that); return *this;} 

    const double* const GetUAnchor() {return u_anchor_;}
    const double* const GetVAnchor() {return v_anchor_;}
    const double* const GetWAnchor() {return w_anchor_;}

    void InitConfigurator(const double r,
                          const double theta,
                          const double phi,
                          const double diameter,
                          const double attach_diameter);
    void Dump() {
      std::cout << std::setprecision(16) << "{" << GetOID() << "," << GetRID() << "," << GetCID() << "}\n"
        << " -> x(" << GetPosition()[0] << ", " << GetPosition()[1] << ", " << GetPosition()[2] << ")\n"
        << " -> u(" << u_anchor_[0] << ", " << u_anchor_[1] << ", " << u_anchor_[2] << ")\n"
        << " -> v(" << v_anchor_[0] << ", " << v_anchor_[1] << ", " << v_anchor_[2] << ")\n"
        << " -> w(" << w_anchor_[0] << ", " << w_anchor_[1] << ", " << w_anchor_[2] << ")\n"
        << "f(" << GetForce()[0] << ", " << GetForce()[1] << ", " << GetForce()[2] << "), "
        << "u(" << GetKineticEnergy() << "), p(" << GetPotentialEnergy() << "), "
        << "d(" << diameter_ << "), attach_diameter(" << attach_diameter_ << ")\n";
    }

    void PrintSPBProperties(int ispb) {
      std::cout << "New Spindle Pole Body: " << GetOID() << std::endl;
      std::cout << std::setprecision(16)
        << "    x(" << GetPosition()[0] << ", " << GetPosition()[1] << ", " << GetPosition()[2] << ")\n"
        << "    u(" << u_anchor_[0] << ", " << u_anchor_[1] << ", " << u_anchor_[2] << ")\n"
        << "    v(" << v_anchor_[0] << ", " << v_anchor_[1] << ", " << v_anchor_[2] << ")\n"
        << "    w(" << w_anchor_[0] << ", " << w_anchor_[1] << ", " << w_anchor_[2] << ")\n"
        << "    translational gamma: " << gamma_tra_ << "\n"
        << "    rotational gamma: " << gamma_rot_ << "\n";
    }

    // Special functions
    void UpdateSPBRefVecs();
    void UpdateSPBDragConstants();

    // Movement along spherical boundary
    void UpdatePositionMP();
};


#include "species.h"
class SpindlePoleBodySpecies : public Species<SpindlePoleBody> {

  protected:

  public:
    SpindlePoleBodySpecies() : Species() {
      SetSID(SID::spb);
    }
    SpindlePoleBodySpecies(int n_members, system_parameters *params, space_struct *space, long seed) : Species(n_members, params, space, seed) {
      SetSID(SID::spb);
    }
    ~SpindlePoleBodySpecies() {}
    SpindlePoleBodySpecies(const SpindlePoleBodySpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }

    void Configurator();
    void ConfiguratorSpindle(int ispb, al_set *anchors);

};

#endif
