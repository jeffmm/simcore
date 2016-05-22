// Contains the main properties structs to run this and control the program flow, etc

#ifndef BUFFMD_PROPERTIES_H_
#define BUFFMD_PROPERTIES_H_

enum ForceType : unsigned char {
    BRUTEFORCE,
    FCELLS,
    FNEIGHBORS_ALLPAIRS,
    FNEIGHBORS_CELL
};

struct _properties {
    _properties(int pCellUpdFreq,
                ForceType pFT) : cell_update_freq_(pCellUpdFreq),
                                 scheme_(pFT) {}
    
    int nparticles_;
    int nspecies_;
    double box_[3];
    double skin_;
    double dt_;
    
    const ForceType scheme_ = FCELLS;
    const int cell_update_freq_ = 4;
};
typedef struct _properties properties_t;

#endif