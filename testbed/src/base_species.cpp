#include "base_species.h"

BaseSpecies::BaseSpecies() {}
BaseSpecies::BaseSpecies(double pRcut,
            int pSid,
            std::string pName) : sid_(pSid), rcut_(pRcut) {
    nparticles_ = 0;
    name_ = pName;
    meff_ = 1.0;
}


BaseSpecies::~BaseSpecies() {

}


double
BaseSpecies::Ukin(std::vector<particle*>* particles) {
    buffmd::unused(particles);
    ukin_ = 0.0;
    return ukin_;
}


void
BaseSpecies::print() {
    printf("{name:%s}, {rcut: %f}, {sid: %d}\n",
            name_.c_str(), rcut_, sid_);
}


void
BaseSpecies::dump() {
    // Dump all information from BaseSpecies
    std::cout << "\tBase Species Dump: \n\t";
    BaseSpecies::print();
}


void
BaseSpecies::checkParticles() {

}


int
BaseSpecies::getSid() {
    return sid_;
}


int
BaseSpecies::getNParticles() {
    return nparticles_;
}


double
BaseSpecies::getRcut() {
    return rcut_;
}


double
BaseSpecies::getMeff() {
    return meff_;
}

