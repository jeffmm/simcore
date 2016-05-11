// Most basic species, needed for object aliasing

#ifndef BUFFMD_BASE_SPECIES_H_
#define BUFFMD_BASE_SPECIES_H_

#include <vector>
#include <string>

#include "helpers.h"
#include "particle.h"

class BaseSpecies {
public:

    BaseSpecies();
    BaseSpecies(double pRcut = 0.0,
                int pSid = -1,
                std::string pName = "");
    virtual ~BaseSpecies();

    virtual void print();
    virtual void dump();

    virtual void checkParticles();

    virtual particle* addParticle(int, int) = 0;
    virtual std::pair<double, double> Ukin(std::vector<particle*>* particles);

    int getSid();
    int getNParticles();
    double getRcut();
    double getMeff();
    
protected:

    int sid_;
    int nparticles_;
    double rcut_;
    double ukin_;
    double temp_;
    double meff_;
    std::string name_;
};

#endif

