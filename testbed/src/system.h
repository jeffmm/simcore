// Base system for others to derive from

#ifndef BUFFMD_SYSTEM_H_
#define BUFFMD_SYSTEM_H_

#include <vector>
#include <unordered_map>
#include <algorithm>

#include "properties.h"
#include "particle.h"
#include "species.h"
#include "cell_list.h"
#include "cell_list_adj.h"
#include "neighbor_list.h"
#include "neighbor_list_cell.h"
#include "potential_manager.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

class SystemArch {
public:
    SystemArch(properties_t* pProperties);
    ~SystemArch();
    
    // Setters (adders)
    void addSpecies(int sid, BaseSpecies* new_species);
    void addPotential(int pid1, int pid2, PotentialBase* newPot);
    int addParticle(int idx);
    
    // Getters
    BaseSpecies* getSpecies(int idx);
    PotentialBase* getPotential(int pid1, int pid2);
    particle* getParticle(int idx);
    int nParticles();
    
    // Shared data structure uses
    void initMP();
    void flattenParticles();
    void generateCellList();
    void updateCellList();
    void generateNeighborList();
    void generateNeighborListCell();
    
    // IO routines and information print
    void output(FILE* erg, FILE* traj, int nfi);
    void dump();
    void dumpPotentials();
    void checkConsistency();
    void statistics(int pNsteps);
    
    // Force calculation and integration
    void forceBrute();
    void forceCellsMP();
    void forceNeighAP();
    void forceNeighCell();
    void velverlet();
    void calcPotential(int psid1, int psid2, double x[3], double y[3], double* fpote);
    std::pair<double, double> ukin();

protected:
    int nparticles_;
    int nsys_;
    int next_pid_;
    int nthreads_;
    int ndim_;

    double box_[3];
    double upot_;
    double ukin_;
    double temperature_;
    double dt_;
    double skin_;
    
    const double kboltz = 0.0019872067;

    properties_t* system_properties_;
    CellList cell_list_;
    //CellListAdj cell_list_adj_;
    
    NeighborList neighbor_list_;
    NeighborListCell neighbor_list_cell_;
    PotentialManager potential_manager_;
    std::unordered_map<int, BaseSpecies*> species_;
    std::vector<particle*> particles_;
    std::vector<double> frc_;
    std::vector<std::pair<int, int>> ranges_;
};

#endif
