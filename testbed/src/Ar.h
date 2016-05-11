// Argon class

#ifndef BUFFMD_AR_H_
#define BUFFMD_AR_H_

#include "particle.h"

class Ar : public particle {
public:
    template<typename... ARGS>
    Ar(ARGS&&... args) : particle{ std::forward<ARGS>(args)... } {}

    particle *getParticle() {
        return this;
    }
};

#endif
