// Species keeps track of a bunch of particles, common properties

#ifndef BUFFMD_SPECIES_H_
#define BUFFMD_SPECIES_H_

#include <vector>
#include <assert.h>

#include "helpers.h"
#include "particle.h"
#include "base_species.h"

template<class T>
class Species : public BaseSpecies {
public:

    typedef std::shared_ptr<T> ptype;
    Species(double pSigma = 0.0,
            double pRcut = 0.0,
            int pSid = -1) : BaseSpecies(pRcut, pSid), sigma_(pSigma) {
    }

    // Set function for Species variables
    virtual void setSpecies(double pSigma,
                            double pRcut,
                            int pSid,
                            std::string pName) {
        sigma_ = pSigma;
        rcut_ = pRcut;
        sid_ = pSid;
        name_ = pName;
    }

    // Print only the variables I have
    virtual void print() {
        printf("{sigma: %f}, {nparticles: %d}\n",
                sigma_, nparticles_);
    }

    // Recursively dump information
    virtual void dump() {
        // Dump all information from Species
        std::cout << "\tSpecies Dump: \n\t";
        Species::print();
        BaseSpecies::dump();
    }

    // Add a particle to particles
    virtual particle* addParticle(int sid, int pid) {
        auto part = particleFactory<T>(pid, sid, name_);
        nparticles_++;
        return part;
    }

protected:

    double sigma_;
};

#endif
