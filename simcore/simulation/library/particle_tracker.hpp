#include "auxiliary.hpp"
#include "cell_list.hpp"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class ParticleTracker {
 private:
  system_parameters* params_;
  int n_dim_, n_per_;
  double cell_length_1d_;
  int n_cells_1d_;
  Cell*** clist_;
  std::vector<ix_pair>* nlist_;
  std::vector<Object*>* objs_;
  void AllocateCellList();
  void DeallocateCellList();

 public:
  ParticleTracker() {}
  void Init(system_parameters* p, std::vector<Object*>* o,
            std::vector<ix_pair>* n);
  void CreatePairsCellList();
  void CreatePartialPairsCellList(std::vector<Object*> ixs, int n_interactors);
  void AddToCellList(std::vector<Object*> ixs, int n_interactors);
  void AssignCells();
  void ClearCells();
  void Clear();
  void CreatePairs();
  double GetCellLength();
};
