// Argon species derived from species

#ifndef BUFFMD_ARSPECIES_H_
#define BUFFMD_ARSPECIES_H_

#include "Ar.h"
#include "species.h"

class ArSpecies : public Species<Ar> {
public:

    template<typename... ARGS>
    ArSpecies(double pMass = 0.0,
              double pEps = 0.0,
              ARGS&&... args) : Species { std::forward<ARGS>(args)... },
                                mass(pMass), eps(pEps) {}

    template<typename... ARGS>
    void setSpecies(double pMass,
                    double pEps,
                    ARGS&&... args) {
        mass = pMass;
        eps = pEps;
        // Call the next initializer along
        Species::setSpecies(std::forward<ARGS>(args)...);
        meff_ = mass * mvsq2e;
    }

    virtual void print() {
        printf("{mass: %f}, {eps: %f}\n", mass, eps);
    }

    virtual double Ukin(std::vector<particle*>* particles) {
        ukin_ = 0.0;
        for (auto it : *particles) {
            if (it->sid != sid_) continue;
            ukin_ += it->v[0]*it->v[0] + it->v[1]*it->v[1]
                     + it->v[2]*it->v[2];
        }
        ukin_ *= 0.5 * mvsq2e * mass;
        return ukin_;
    }

    virtual void dump() {
        printf("%s\n\t", name_.c_str());
        ArSpecies::print();
        Species::dump();
    }

    virtual particle* addParticle(int sid, int pid) {
        return Species::addParticle(sid, pid);
    }

    // Constants for all argon particles
    double mass;
    double eps;
    const double mvsq2e = 2390.05736153349;
};

#endif
