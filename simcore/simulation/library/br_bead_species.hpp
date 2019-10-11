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
  BrBeadSpecies() : Species() {}
  void Init(system_parameters *params, species_base_parameters *sparams,
            space_struct *space);
  void UpdatePositions();
};

#endif
