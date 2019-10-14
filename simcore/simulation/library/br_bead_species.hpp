#ifndef _SIMCORE_BR_BEAD_SPECIES_H_
#define _SIMCORE_BR_BEAD_SPECIES_H_

#include "species.hpp"
#include "br_bead.hpp"

typedef std::vector<BrBead>::iterator br_bead_iterator;
typedef std::vector<
    std::pair<std::vector<BrBead>::iterator, std::vector<BrBead>::iterator>>
    br_bead_chunk_vector;

class BrBeadSpecies : public Species<BrBead, species_id::br_bead> {
public:
  BrBeadSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void UpdatePositions();
};

#endif
