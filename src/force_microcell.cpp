// Implementation for force microcells

#include "force_microcell.h"

void
ForceMicrocell::InitMP() {
    printf("ForceMicrocell::InitMP\n");
    microcell_list_.Init(space_, skin_);
    //microcell_list_.LoadFlatSimples(simples_);
    //microcell_list_.CreateSubstructure(1.0);
    //XXX stuff
}


void
ForceMicrocell::Finalize() {
    //XXX
    initialized_ = false;
}


void
ForceMicrocell::UpdateScheme() {
    // XXX stuff
}


void
ForceMicrocell::Interact() {
    // Stuff!
}
