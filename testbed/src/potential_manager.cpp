// Potential Manager implementation

#include <algorithm>

#include "potential_manager.h"

PotentialManager::PotentialManager(void) {

}

PotentialManager::~PotentialManager() {
    for (auto pot : potentials) {
        delete(pot.second);
    }
}

void
PotentialManager::print() {
    for(auto& x : potentials) {
        printf("{%d,%d} ->\n", x.first.first, x.first.second);
        x.second->print();
    }
}

void
PotentialManager::addPotential(int ptid1, int ptid2, PotentialBase* newPot) {
    auto key = std::make_pair(ptid1, ptid2);
    potentials[key] = newPot;
}

PotentialBase*
PotentialManager::getPotential(int ptid1, int ptid2) {
    auto key = std::make_pair(ptid1, ptid2);
    return potentials[key];
}


void
PotentialManager::calculatePotential(int psid1, int psid2, double* x, double* y, double* fpote) {
    auto key = std::make_pair(psid1, psid2);
    auto pot = potentials.find(key);
    if(pot != potentials.end()) {
        potentials[key]->CalcPotential(x, y, fpote);
    } else {
        std::fill(fpote, fpote + 4, 0.0);
    }
}
