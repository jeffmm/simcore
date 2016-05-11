// Magical potential manager for the system

#ifndef BUFFMD_POTENTIAL_MANAGER_H_
#define BUFFMD_POTENTIAL_MANAGER_H_

#include <map>
#include <memory>

#include "helpers.h"
#include "potential_base.h"
#include "lennard_jones_12_6.h"

class PotentialManager {
public:
    PotentialManager();
    ~PotentialManager();
    void print();
    void addPotential(int ptid1, int ptid2, PotentialBase* newPot);
    PotentialBase* getPotential(int ptid1, int ptid2);
    void calculatePotential(int psid1, int psid2, double* x, double* y, double* fpote);

protected:
    std::map<std::pair<int, int>, PotentialBase*> potentials;
};

#endif
