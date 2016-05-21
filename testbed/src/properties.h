// Contains the main properties structs to run this and control the program flow, etc

#ifndef BUFFMD_PROPERTIES_H_
#define BUFFMD_PROPERTIES_H_

struct _properties {
    _properties(int pCellUpdFreq,
                bool pUseCells) : cell_update_freq_(pCellUpdFreq),
                                  use_cells_(pUseCells) {}
    
    int nparticles_;
    int nspecies_;
    double box_[3];
    double skin_;
    double dt_;
    
    const bool use_cells_ = true;
    const int cell_update_freq_ = 4;
};
typedef struct _properties properties_t;

#endif