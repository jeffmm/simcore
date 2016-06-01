// Contains the main properties structs to run this and control the program flow, etc

#ifndef BUFFMD_PROPERTIES_H_
#define BUFFMD_PROPERTIES_H_

enum ForceType : unsigned char {
    BRUTEFORCE,
    FCELLS,
    FCELLSADJ,
    FCELLSADJ_INT,
    FNEIGHBORS_ALLPAIRS,
    FNEIGHBORS_CELL
};

struct _properties {
    _properties(int pCellUpdFreq,
                ForceType pFT) : cell_update_freq_(pCellUpdFreq),
                                 scheme_(pFT) {}
    
    int nparticles_ = 0;
    int nspecies_ = 0;
    int ndim_ = 3;
    double box_[3] = {0.0, 0.0, 0.0};;
    double skin_ = 0.0;
    double dt_ = 0.0;
    
    const int cell_update_freq_ = 4;
    const ForceType scheme_ = FCELLS;
};
typedef struct _properties properties_t;

#endif
