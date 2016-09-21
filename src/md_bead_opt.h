#ifndef _SIMCORE_MD_BEAD_OPT_H_
#define _SIMCORE_MD_BEAD_OPT_H_

#include "md_bead.h"
#include "output_manager.h"

class MDBeadOpt : public MDBead {
  public:
    MDBeadOpt(system_parameters *params, space_struct *space, 
        long seed, SID sid) : MDBead(params, space, seed, sid) {}
    ~MDBeadOpt() {}

  //void WritePosit(std::ofstream &op);
  //void ReadPosit(std::ifstream &ip);

};

class MDBeadOptSpecies : public Species<MDBeadOpt> {
  private:

  protected:

  public:
    MDBeadOptSpecies() : Species() {
      SetSID(SID::md_bead_opt);
    }

    MDBeadOptSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      SetSID(SID::md_bead_opt);
    }

    ~MDBeadOptSpecies() {}
    MDBeadOptSpecies(const MDBeadOptSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) { //FIXME not sure if this is right
      SpeciesBase::operator=(that);
      return *this;
    }
    void Init() {
      Species::Init();
    }

    // Configurations
    void Configurator(); 

    //void WritePosits();
    //void ReadPosits();
};


#endif //_SIMCORE_MD_BEAD_OPT_H_
