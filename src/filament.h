#ifndef _SIMCORE_FILAMENT_H_
#define _SIMCORE_FILAMENT_H_

#include "site.h"
#include "bond.h"
#include "species.h"
#include "auxiliary.h"
#include "wca.h"

class Filament : public Composite<Site,Bond> {

  private:
    int n_bonds_,
        n_sites_,
        dynamic_instability_flag_,
        force_induced_catastrophe_flag_;
    double max_length_,
           min_length_,
           max_child_length_,
           child_length_,
           gamma_par_,
           gamma_perp_,
           rand_sigma_par_,
           rand_sigma_perp_,
           v_poly_,
           v_depoly_,
           p_s2g_,
           p_s2p_,
           p_p2s_,
           p_p2g_,
           p_g2s_,
           p_g2p_,
           tip_force_,
           **gamma_inverse_,
           *p_vec_, //n_sites-1
           *g_mat_lower_, //n_sites-2
           *g_mat_upper_, //n_sites-2
           *g_mat_diag_, //n_sites-1
           *det_t_mat_, //n_sites+1
           *det_b_mat_, //n_sites+1
           *g_mat_inverse_, //n_sites-2
           *k_eff_, //n_sites-2
           *h_mat_diag_, //n_sites-1
           *h_mat_upper_, //n_sites-2
           *h_mat_lower; //n_sites-2
    poly_state_t poly_state_;
    void UpdateSiteBondPositions();
    void SetDiffusion();
    void UpdateOrientation();
    void GetBodyFrame();
    void AddRandomDisplacement();
    void ApplyForcesTorques();
    void DynamicInstability();

  public:
    Filament(system_parameters *params, space_struct * space, long seed, SID sid) 
      : Composite(params, space, seed, sid) {
        //bond_ptr_ = new Bond(params, space, gsl_rng_get(rng_.r), GetSID()) ;
        length_ = params->rod_length;
        diameter_ = params->rod_diameter;
        max_length_ = params->max_rod_length;
        min_length_ = params->min_rod_length;
        max_child_length_ = 0.5*params->cell_length;
        dynamic_instability_flag_ = params->dynamic_instability_flag;
        force_induced_catastrophe_flag_ = params->force_induced_catastrophe_flag;
        p_g2s_ = params->f_grow_to_shrink*delta_;
        p_g2p_ = params->f_grow_to_pause*delta_;
        p_s2p_ = params->f_shrink_to_pause*delta_;
        p_s2g_ = params->f_shrink_to_grow*delta_;
        p_p2s_ = params->f_pause_to_shrink*delta_;
        p_p2g_ = params->f_pause_to_grow*delta_;
        v_depoly_ = params->v_depoly;
        v_poly_ = params->v_poly;

        n_bonds_ = (int) ceil(length_/max_child_length_);
        n_sites_ = n_bonds_+1;
        child_length_ = length_/n_bonds_;
        // Initialize sites
        for (int i=0; i<n_sites_; ++i) {
          Site s(params, space, gsl_rng_get(rng_.r), GetSID());
          s.SetCID(GetCID());
          elements_.push_back(s);
        }
        // Initialize bonds
        for (int i=0; i<n_bonds_; ++i) {
          Bond b(params, space, gsl_rng_get(rng_.r), GetSID());
          b.SetCID(GetCID());
          b.SetRID(GetRID());
          //b.InitOID();
          v_elements_.push_back(b);
        }
        //Allocate control structures
        //TODO: JMM Consider another method of storing these arrays,
        //as it is a lot of memory for each filament to store. This
        //may be the only way to do this if we want to use parallel
        int max_sites = (int) ceil(max_length_/min_child_length_) + 1; 
        p_vec_ = new double[max_sites-1];
        g_mat_lower_ = new double[max_sites-2];
        g_mat_upper_ = new double[max_sites-2];
        g_mat_diag_ = new double[max_sites-1];
        det_t_mat_ = new double[max_sites+1];
        det_b_mat_ = new double[max_sites+1];
        g_mat_inverse_ = new double[max_sites-2];
        k_eff_ = new double[max_sites-2];
        h_mat_diag_ = new double [max_sites-1];
        h_mat_upper_ = new double [max_sites-2];
        h_mat_lower_ = new double [max_sites-2];
        gamma_inverse_ = new double*[max_sites];
        for(int i=0; i<max_sites; ++i)
          gamma_inverse_[i] = new double[n_dim_*n_dim_];
      }
    ~Filament() {
      //Clean up control structures
      delete[] p_vec_;
      delete[] g_mat_lower_;
      delete[] g_mat_upper_;
      delete[] g_mat_diag_;
      delete[] det_t_mat_;
      delete[] det_b_mat_;
      delete[] g_mat_inverse_;
      delete[] k_eff_;
      delete[] h_mat_diag_;
      delete[] h_mat_upper_;
      delete[] h_mat_lower_;
      int max_sites = (int) ceil(max_length_/min_child_length_)+1;
      for (int i=0; i<max_sites; ++i)
        delete[] gamma_inverse_[i];
      delete[] gamma_inverse_;
    }
    Filament(const Filament& that) : Composite(that) {
      n_bonds_=that.n_bonds_;
      n_sites_ = that.n_sites_;
      dynamic_instability_flag_ = that.dynamic_instability_flag_;
      max_length_ = that.max_length_;
      min_length_ = that.min_length_;
      max_child_length_ = that.max_child_length_;
      gamma_par_ = that.gamma_par_;
      gamma_perp_ = that.gamma_perp_;
      gamma_rot_ = that.gamma_rot_;
      rand_sigma_par_ = that.rand_sigma_par_;
      rand_sigma_perp_ = that.rand_sigma_perp_;
      rand_sigma_rot_ = that.rand_sigma_rot_;
      v_poly_ = that.v_poly_;
      v_depoly_ = that.v_depoly_;
      p_s2g_ = that.p_s2g_;
      p_s2p_ = that.p_s2p_;
      p_p2s_ = that.p_p2s_;
      p_p2g_ = that.p_p2g_;
      p_g2s_ = that.p_g2s_;
      p_g2p_ = that.p_g2p_;
      tip_force_ = that.tip_force_;
      std::copy(that.body_frame_, that.body_frame_+6, body_frame_);
      gamma_inverse_ = that.gamma_inverse_;
      p_vec_ = that.p_vec_;
      g_mat_lower_ = that.g_mat_lower_;
      g_mat_upper_ = that.g_mat_upper_;
      g_mat_diag_ = that.g_mat_diag_;
      h_mat_diag_ = that.h_mat_diag_;
      h_mat_lower_ = that.h_mat_lower_;
      h_mat_upper_ = that.h_mat_upper_;
      k_eff_ = that.k_eff_;
      g_mat_inverse_ = that.g_mat_inverse_;
      det_t_mat_ = that.det_t_mat_;
      det_b_mat_ = that.det_b_mat_;
      poly_state_ = that.poly_state_;
    }
    Filament& operator=(Filament const& that) {
      Composite::operator=(that); 
      n_bonds_=that.n_bonds_;
      dynamic_instability_flag_ = that.dynamic_instability_flag_;
      max_length_ = that.max_length_;
      min_length_ = that.min_length_;
      max_child_length_ = that.max_child_length_;
      gamma_par_ = that.gamma_par_;
      gamma_perp_ = that.gamma_perp_;
      gamma_rot_ = that.gamma_rot_;
      rand_sigma_par_ = that.rand_sigma_par_;
      rand_sigma_perp_ = that.rand_sigma_perp_;
      rand_sigma_rot_ = that.rand_sigma_rot_;
      v_poly_ = that.v_poly_;
      v_depoly_ = that.v_depoly_;
      p_s2g_ = that.p_s2g_;
      p_s2p_ = that.p_s2p_;
      p_p2s_ = that.p_p2s_;
      p_p2g_ = that.p_p2g_;
      p_g2s_ = that.p_g2s_;
      p_g2p_ = that.p_g2p_;
      tip_force_ = that.tip_force_;
      std::copy(that.body_frame_, that.body_frame_+6, body_frame_);
      poly_state_ = that.poly_state_;
      return *this;
    } 
    virtual void Init();
    virtual void Integrate();
    virtual double const * const GetDrTot();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePositionMP();
};

class FilamentSpecies : public Species<Filament> {
  protected:
    //void InitPotentials(system_parameters *params);
  public:
    FilamentSpecies(int n_members, system_parameters *params, 
        space_struct *space, long seed) 
      : Species(n_members, params, space, seed) {
      SetSID(SID::br_rod);
      //InitPotentials(params);
    }
    ~FilamentSpecies() {}
    FilamentSpecies(const FilamentSpecies& that) : Species(that) {}
    Species& operator=(Species const& that) {
      SpeciesBase::operator=(that);
      return *this;
    }
    void Init() {
      Species::Init();
    }

};

#endif // _SIMCORE_FILAMENT_H_
