// Implementation for force microcells

#include "force_microcell.h"

void
ForceMicrocell::InitMP() {
    microcell_list_.Init(space_, skin_);
    microcell_list_.LoadFlatSimples(simples_);
    //XXX stuff
}


void
ForceMicrocell::UpdateScheme() {
    // XXX stuff
}


void
ForceMicrocell::Interact() {
    // Stuff!
}
