// Base system for others to derive from

#ifndef BUFFMD_SYSTEM_H_
#define BUFFMD_SYSTEM_H_

#include <vector>
#include <unordered_map>
#include <algorithm>

#include "particle.h"
#include "species.h"
#include "cell_list.h"
#include "potential_manager.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

class SystemArch {
public:
    SystemArch(int nparticles, double pBox, double pDt);
    ~SystemArch();
    
    void addSpecies(int sid, BaseSpecies* new_species);
    void addPotential(int pid1, int pid2, PotentialBase* newPot);
    int addParticle(int idx);
    
    BaseSpecies* getSpecies(int idx);
    PotentialBase* getPotential(int pid1, int pid2);
    particle* getParticle(int idx);
    
    void flattenParticles();
    void generateCellList();
    void updateCellList();
    void dump();
    void dumpPotentials();
    void checkConsistency();
    void forceMP();
    void velverlet();
    void output(FILE* erg, FILE* traj, int nfi);
    void calcPotential(int psid1, int psid2, double x[3], double y[3], double* fpote);

    void testCellLists();

    std::pair<double, double> ukin();
    int nParticles();

protected:
    int nparticles_;
    int nsys_;
    int next_pid_;
    int nthreads_;

    double box_;
    double upot_;
    double ukin_;
    double temperature_;
    double dt_;

    CellList cell_list_;
    PotentialManager potential_manager_;
    std::unordered_map<int, BaseSpecies*> species_;
    std::vector<particle*> particles_;
    std::vector<double> frc_;
    std::vector<std::pair<int, int>> ranges_;
};

#endif
