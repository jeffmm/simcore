// Main neighbor list implementation, should be derived from

// NOTE: Implements all pairs search by default!!!!

#ifndef BUFFMD_NEIGHBOR_LIST_H_
#define BUFFMD_NEIGHBOR_LIST_H_

#include <vector>

#include "particle.h"

struct _neighbor {
    int idx_;
    double value_;
};
typedef struct _neighbor neighbor_t;

typedef std::vector<neighbor_t> nl_list;

class NeighborList {
public:
    
    NeighborList() {}
    virtual ~NeighborList();
   
    virtual void print();
    virtual void CreateNeighborList(int pN, double pRcut, double pSkin, double pBox[3]);
    virtual void CheckNeighborList(std::vector<particle*>* particles);
    virtual void UpdateNeighborList(std::vector<particle*>* particles);
    
    void AllPairsUpdate(std::vector<particle*>* particles);
    
    virtual unsigned long GetMemoryFootprint();
    
    // Getters
    int GetNUpdates() { return n_updates_; }
    const nl_list* GetNeighbors() {
        return neighbors_;
    }
    
    
protected:
    
    // Inputs
    int nparticles_;
    double rcut_;
    double skin_;
    double box_[3];
    
    // Computed quantities
    double rcut2_; // Cutoff radius squared
    double skin2_; // Skin radius squared
    double rcs2_; // (cut + skin)^2
    double half_skin2_; // Half the skin distance squared (needed for updates)
    double boxby2_[3];
    
    // Internal flags, statistics
    bool nl_update_;
    int n_updates_;
    
    // The neighbor list itself, hooray!
    // Might be inefficient in terms of access time, etc, because of dynamic array
    nl_list* neighbors_;
};

#endif
