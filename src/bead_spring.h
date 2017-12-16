#ifndef _SIMCORE_BEAD_SPRING_H_
#define _SIMCORE_BEAD_SPRING_H_

#include "species.h"
#include "mesh.h"

#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

class BeadSpring : public Mesh {

  private:
    int stoch_flag_,
        eq_steps_,
        eq_steps_count_ = 0;
    double max_bond_length_,
           bond_spring_,
           persistence_length_,
           rand_sigma_,
           driving_factor_;
    std::vector<double> cos_thetas_;
    void SetDiffusion();
    void GenerateProbableOrientation();
    void CalculateAngles();
    void CalculateTangents();
    void AddRandomForces();
    void AddBendingForces();
    void AddSpringForces();
    void UpdateSitePositions();
    void UpdateSiteOrientations();
    void ApplyForcesTorques();
    void ApplyInteractionForces();
    void SetParameters();
    void InitElements();
    void InsertBeadSpring(bool force_overlap=false);
    void InsertFirstBond();
    void UpdateAvgPosition();
    void ReportAll();

  public:
    BeadSpring();
    virtual void Init(bool force_overlap = false);
    virtual void InsertAt(double *pos, double *u);
    virtual bool CheckBounds(double buffer = 0);
    //void DiffusionValidationInit();
    virtual void Integrate();
    virtual void Draw(std::vector<graph_struct*> * graph_array);
    virtual void UpdatePosition();
    double const GetLength() { return length_;}
    double const GetDriving() {return driving_factor_;}
    double const GetPersistenceLength() {return persistence_length_;}
    double const GetBondLength() {return bond_length_;}
    int const GetNBonds() {return n_bonds_;}
    void GetAvgOrientation(double * au);
    void GetAvgPosition(double * ap);
    std::vector<double> const * const GetThetas() {
      return &cos_thetas_;
    }
    double GetTipZ() {
      return sites_[n_sites_-1].GetOrientation()[n_dim_-1];
    }
    double const * const GetHeadPosition() {
      return sites_[n_sites_-1].GetPosition();
    }
    double const * const GetTailPosition() {
      return sites_[0].GetPosition();
    }
    double const * const GetTailOrientation() {
      return sites_[0].GetOrientation();
    }
    void AddTorqueTail(double * t) {
      bonds_[0].AddTorque(t);
    }
    void AddForceTail(double *f) {
      sites_[0].AddForce(f);
    }
    void WritePosit(std::fstream &oposit);
    void ReadPosit(std::fstream &iposit);
    void WriteSpec(std::fstream &ospec);
    void ReadSpec(std::fstream &ispec);
    void WriteCheckpoint(std::fstream &ocheck);
    void ReadCheckpoint(std::fstream &icheck);
    void ScalePosition();
    double const GetVolume();
};

typedef std::vector<BeadSpring>::iterator bs_iterator;
typedef std::vector<std::pair<std::vector<BeadSpring>::iterator, std::vector<BeadSpring>::iterator> > bead_spring_chunk_vector;

class BeadSpringSpecies : public Species<BeadSpring> {
  protected:
    // Analysis structures
    double e_bend_,
           tot_angle_,
           mse2e_,
           mse2e2_;
    int **theta_histogram_;
    int time_,
        n_bins_,
        n_samples_;
    std::fstream spiral_file_,
                 theta_file_,
                 mse2e_file_;
  public:
    BeadSpringSpecies() : Species() {
      SetSID(species_id::bead_spring);
    }
    void Init(system_parameters *params, space_struct *space, long seed) {
      Species::Init(params, space, seed);
      sparams_ = &(params_->bead_spring);
      if (params_->bead_spring.packing_fraction>0) {
        if (params_->bead_spring.length <= 0) {
          error_exit("Packing fraction with polydisperse lengths not implemented yet\n");
        }
        if (params_->n_dim == 2) {
          double fil_vol = params_->bead_spring.length*params_->bead_spring.diameter+0.25*M_PI*SQR(params_->bead_spring.diameter);
          sparams_->num = params_->bead_spring.packing_fraction*space_->volume/fil_vol;
        }
        else {
          double fil_vol = 0.25*M_PI*SQR(params_->bead_spring.diameter)*params_->bead_spring.length+M_PI*CUBE(params_->bead_spring.diameter)/6.0;
          sparams_->num = params_->bead_spring.packing_fraction*space_->volume/fil_vol;
        }
      }
    }
    void InitAnalysis();
    void InitThetaAnalysis();
    void InitMse2eAnalysis();
    void RunAnalysis();
    void RunThetaAnalysis();
    void RunMse2eAnalysis();
    void FinalizeAnalysis();
    void FinalizeMse2eAnalysis();
    void FinalizeThetaAnalysis();
    void UpdatePositions() {
#ifdef ENABLE_OPENMP
      int max_threads = omp_get_max_threads();
      bead_spring_chunk_vector chunks;
      chunks.reserve(max_threads); 
      size_t chunk_size= members_.size() / max_threads;
      bs_iterator cur_iter = members_.begin();
      for(int i = 0; i < max_threads - 1; ++i) {
        bs_iterator last_iter = cur_iter;
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
      for (bs_iterator it=members_.begin(); it!=members_.end(); ++it) 
        it->UpdatePosition();
#endif

    }
};

#endif // _SIMCORE_BEAD_SPRING_H_
