#include <simcore/br_bead_species.hpp>

BrBeadSpecies::BrBeadSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::br_bead);
}
void BrBeadSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
}

void BrBeadSpecies::UpdatePositions() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  br_bead_chunk_vector chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = members_.size() / max_threads;
  br_bead_iterator cur_iter = members_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    br_bead_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i)
      for (auto it = chunks[i].first; it != chunks[i].second; ++it)
        it->UpdatePosition();
  }
#else
  for (br_bead_iterator it = members_.begin(); it != members_.end(); ++it)
    it->UpdatePosition();
#endif
}
