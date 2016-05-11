#ifndef BUFFMD_PARTICLE_H_
#define BUFFMD_PARTICLE_H_

// Need to include information about namespace std, stupid thing
#include <cstring>
#include <iostream>
#include <memory>

#include "helpers.h"

// Base Particle class
class particle {
public:
    particle(const int &pPid = -1,
             const int &pType = -1,
             std::string pName = "") : pid(pPid), sid(pType) {
        name = pName;
        init();
    }

    virtual ~particle() {

    }

    void init() {
        memset(x, 0, sizeof(x));
        memset(v, 0, sizeof(v));
        memset(f, 0, sizeof(f));
    }
    void print() {
        printf("pid: %d, sid: %d\n", pid, sid);
        printf("x: {%20.8f, %20.8f, %20.8f}\n", x[0], x[1], x[2]);
        printf("v: {%20.8f, %20.8f, %20.8f}\n", v[0], v[1], v[2]);
        printf("f: {%20.8f, %20.8f, %20.8f}\n", f[0], f[1], f[2]);
    }
    void setXYZ(double rx, double ry, double rz) {
        x[0] = rx;
        x[1] = ry;
        x[2] = rz;
    }
    void setV(double vx, double vy, double vz) {
        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    }

    void zeroF() {
        memset(f, 0, sizeof(f));
    }

    virtual particle *getParticle() = 0;
    
    bool operator < (const particle& part) const {
        return (sid < part.sid);
    }

    // For now, just use an array of doubles
    double x[3]; // position
    double v[3]; // velocity
    double f[3]; // force

    const int pid; // particle id
    const int sid; //species id

    std::string name;
};


// Particle Factory using templates
template<typename T, typename...ARGS, typename = typename std::enable_if<std::is_base_of<particle, T>::value>::type>
T* particleFactory(ARGS&&... args) {
    T* part{ new T{ std::forward<ARGS>(args)...} };

    return part;
}

#endif
