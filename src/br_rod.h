#ifndef _SIMCORE_BR_ROD_H_
#define _SIMCORE_BR_ROD_H_

#include "anchor_list_generic.h"
#include "site.h"
#include "bond.h"
#include "species.h"
#include "auxiliary.h"
#include "wca.h"

#include <unordered_map>

class BrRod : public Composite<Site,Bond> {

  private:
    bool rod_diffusion_;
    bool rod_fixed_;
    int n_bonds_,
        stabilization_state_,
        diffusion_validation_flag_;
    double max_length_,
           min_length_,
           max_child_length_,
           child_length_,
           gamma_par_,
           gamma_perp_,
           gamma_rot_,
           rand_sigma_par_,
           rand_sigma_perp_,
           rand_sigma_rot_,
           f_stabilize_fr_,
           f_stabilize_fc_,
           f_stabilize_vg_,
           f_stabilize_vs_,
           body_frame_[6];
    poly_state_t poly_state_;
    void UpdateSitePositions();
    void UpdateBondPositions();
    void UpdateAnchors();
    void SetDiffusion();
    void UpdateOrientation();
    void GetBodyFrame();
    void AddRandomDisplacement();
    void ApplyForcesTorques();

  public:
    BrRod(system_parameters *params, space_struct * space, long seed, SID sid) 
      : Composite(params, space, seed, sid) {
        length_ = params->rod_length;
        diameter_ = params->rod_diameter;
        max_length_ = params->max_rod_length;
        min_length_ = params->min_rod_length;
        max_child_length_ = 0.5*params->cell_length;
        rod_diffusion_ = params->rod_diffusion == 1 ? true : false;
        rod_fixed_ = params->rod_fixed == 1 ? true : false;
        diffusion_validation_flag_ = params->diffusion_validation_flag;
        // Initialize end sites
        for (int i=0; i<2; ++i) {
          Site s(params, space, gsl_rng_get(rng_.r), GetSID());
          s.SetCID(GetCID());
          elements_.push_back(s);
        }
        // Initialize bonds
        n_bonds_ = (int) ceil(length_/max_child_length_);
        child_length_ = length_/n_bonds_;
        for (int i=0; i<n_bonds_; ++i) {
          Bond b(params, space, gsl_rng_get(rng_.r), GetSID());
          b.SetCID(GetCID());
          b.SetRID(GetRID());
          v_elements_.push_back(b);
        }
      }
    ~BrRod() {}
    BrRod(const BrRod& that) : Composite(that) {
      n_bonds_=that.n_bonds_;
      max_length_ = that.max_length_;
      min_length_ = that.min_length_;
      max_child_length_ = that.max_child_length_;
      gamma_par_ = that.gamma_par_;
      gamma_perp_ = that.gamma_perp_;
      gamma_rot_ = that.gamma_rot_;
      rand_sigma_par_ = that.rand_sigma_par_;
      rand_sigma_perp_ = that.rand_sigma_perp_;
      rand_sigma_rot_ = that.rand_sigma_rot_;
      std::copy(that.body_frame_, that.body_frame_+6, body_frame_);
      poly_state_ = that.poly_state_;
      stabilization_state_ = that.stabilization_state_;
      f_stabilize_fr_ = that.f_stabilize_fr_;
      f_stabilize_fc_ = that.f_stabilize_fc_;
      f_stabilize_vg_ = that.f_stabilize_vg_;
      f_stabilize_vs_ = that.f_stabilize_vs_;
      rod_diffusion_=that.rod_diffusion_;
      rod_fixed_=that.rod_fixed_;
      diffusion_validation_flag_ = that.diffusion_validation_flag_;
    }
    BrRod& operator=(BrRod const& that) {
      Composite::operator=(that); 
      n_bonds_=that.n_bonds_;
      max_length_ = that.max_length_;
      min_length_ = that.min_length_;
      max_child_length_ = that.max_child_length_;
      gamma_par_ = that.gamma_par_;
      gamma_perp_ = that.gamma_perp_;
      gamma_rot_ = that.gamma_rot_;
      rand_sigma_par_ = that.rand_sigma_par_;
      rand_sigma_perp_ = that.rand_sigma_perp_;
      rand_sigma_rot_ = that.rand_sigma_rot_;
      std::copy(that.body_frame_, that.body_frame_+6, body_frame_);
      poly_state_ = that.poly_state_;
      stabilization_state_ = that.stabilization_state_;
      f_stabilize_fr_ = that.f_stabilize_fr_;
      f_stabilize_fc_ = that.f_stabilize_fc_;
      f_stabilize_vg_ = that.f_stabilize_vg_;
      f_stabilize_vs_ = that.f_stabilize_vs_;
      rod_diffusion_ = that.rod_diffusion_;
      rod_fixed_=that.rod_fixed_;
      diffusion_validation_flag_ = that.diffusion_validation_flag_;
      return *this;
    } 
    virtual void Init();
    virtual void Integrate();
    //virtual double const * const GetDrTot();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePositionMP();
    virtual void Dump();

    // Specific functions for configurations
    void InitConfigurator(const double* const x, const double* const u, const double l);
    void InitOriented(const double* const u);

    //void WritePosit(std::fstream &op);
    //void ReadPosit(std::fstream &ip);

    // KMC information
    poly_state_t GetPolyState() {return poly_state_;}
    void SetPolyState(poly_state_t poly_state) {poly_state_=poly_state;}
    int GetStabilizationState(double *f_stabilize_fr,
                              double *f_stabilize_fc,
                              double *f_stabilize_vg,
                              double *f_stabilize_vs) {
      (*f_stabilize_fr) = f_stabilize_fr_;
      (*f_stabilize_fc) = f_stabilize_fc_;
      (*f_stabilize_vg) = f_stabilize_vg_;
      (*f_stabilize_vs) = f_stabilize_vs_;
      return stabilization_state_;
    }
    void SetStabilizationState(const int state,
                               const double f_stab_fr,
                               const double f_stab_fc,
                               const double f_stab_vg,
                               const double f_stab_vs) {
      f_stabilize_fr_ = f_stab_fr;
      f_stabilize_fc_ = f_stab_fc;
      f_stabilize_vg_ = f_stab_vg;
      f_stabilize_vs_ = f_stab_vs;
      stabilization_state_=state;
    }
    double UpdateRodLength(const double delta_length);

};

class BrRodSpecies : public Species<BrRod> {
  protected:
    //void InitPotentials(system_parameters *params);
    bool diffusion_validation_;
    double max_length_;
    double ***orientations_;
    double min_length_;
    int n_dim_,
        nbins_,
        ibin_,
        nvalidate_,
        ivalidate_;
    void ValidateDiffusion();
    void WriteDiffusionValidation(std::string run_name);
  public:
    BrRodSpecies() : Species() {
      SetSID(SID::br_rod);
    }
    BrRodSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      SetSID(SID::br_rod);
      //InitPotentials(params);
      max_length_ = params_->max_rod_length;
    }
    ~BrRodSpecies() {
      if (diffusion_validation_) {
        for (int i=0; i<n_members_; ++i) {
          for (int j=0; j<2; ++j) {
            delete[] orientations_[i][j];
          }
          delete[] orientations_[i];
        }
        delete[] orientations_;
      }
    }
    BrRodSpecies(const BrRodSpecies& that) : Species(that) {
      max_length_ = that.max_length_;
      min_length_ = that.min_length_;
      nbins_ = that.nbins_;
      ibin_ = that.ibin_;
      nvalidate_ = that.nvalidate_;
      ivalidate_ = that.ivalidate_;
      orientations_ = that.orientations_;
    }
    BrRodSpecies& operator=(BrRodSpecies const& that) {
      Species::operator=(that);
      max_length_ = that.max_length_;
      min_length_ = that.min_length_;
      nbins_ = that.nbins_;
      ibin_ = that.ibin_;
      nvalidate_ = that.nvalidate_;
      ivalidate_ = that.ivalidate_;
      orientations_ = that.orientations_;
      return *this;
    }
    //virtual void InitConfig(system_parameters *params, space_struct *space, long seed) {
      //Species::InitConfig(params, space, seed);
    //}
    void Init() {
      Species::Init();
    }
    double const GetMaxLength() {return max_length_;}
    double const GetMinLength() {return min_length_;}

    void UpdatePositionsMP() {
      if (diffusion_validation_ && ivalidate_%nvalidate_ == 0)
        ValidateDiffusion();
      for (auto it=members_.begin(); it != members_.end(); ++it) {
        (*it)->UpdatePositionMP();
      }
      ivalidate_++;
    }
    void WriteOutputs(std::string run_name);

    // Special insertion routine
    void Configurator();
    void ConfiguratorSpindle(int ispb, int spb_oid,
                             const double* const r_spb,
                             const double* const u_spb,
                             const double* const v_spb,
                             const double* const w_spb,
                             al_set *anchors);
    void InsertRandomMT();
    void InsertRTPMT();

    static void CreateTestRod(BrRod **rod,
                              int ndim,
                              std::vector<Simple*>* simples,
                              std::unordered_map<int, int>* oid_position_map,
                              const std::string &filename,
                              const std::string &modulename,
                              const std::string &unitname,
                              const std::string &rodname,
                              int itest);

    static void CreateTestRod(BrRod **rod,
                              int ndim,
                              std::vector<Simple*>* simples,
                              std::unordered_map<int, int>* oid_position_map,
                              YAML::Node *subnode);
                              
};

//TODO Use these in the configurator and not strings
enum class ins_type : unsigned char {
  isotropic,
  polar,
  crystal,
  smectic,
  neumatic
};

#endif // _SIMCORE_BR_ROD_H_
