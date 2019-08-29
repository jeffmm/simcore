#ifndef _SIMCORE_CELL_LIST_H_
#define _SIMCORE_CELL_LIST_H_

#include "auxiliary.hpp"
#include "interaction.hpp"
#include "object.hpp"
//#include <unordered_map>

typedef std::array<int, 3> cell_index;
typedef std::pair<Object*, Object*> interactor_pair;
typedef std::pair<interactor_pair, Interaction> pair_interaction;

class Cell {
 private:
  cell_index index_;
  std::vector<Object*> interactors_;
  std::vector<Cell*> neighbors_;
  std::vector<interactor_pair> interactions_;

 public:
  Cell() {}
  Cell(int i, int j, int k) { index_ = {i, j, k}; }
  void AddInteractor(Object* ix) { interactors_.push_back(ix); }
  void AddNeighbor(Cell* neighbor) { neighbors_.push_back(neighbor); }
  std::vector<Object*>::iterator Begin() { return interactors_.begin(); }
  std::vector<Object*>::iterator End() { return interactors_.end(); }
  int const Count() { return interactors_.size(); }
  Cell* GetCellPtr() { return this; }
  int const GetNInteractions();
  std::vector<interactor_pair> PairInteractions();
};

// typedef std::unordered_map<cell_index, Cell> cell_map;
// typedef std::unordered_map<cell_index, Cell>::iterator cell_map_it;
typedef std::map<cell_index, Cell> cell_map;
typedef std::map<cell_index, Cell>::iterator cell_map_it;

class CellList {
 private:
  int n_dim_;
  int n_periodic_;
  int n_cells_1d_;
  double cell_length_;
  double scaled_cell_length_;
  cell_map cells_;

  void SetNeighbors(int i, int j, int k);
  void InitCells();

 public:
  CellList() {}
  void Init(int n_dim, int n_periodic, double cell_length,
            double system_radius);
  void ClearCells() { cells_.clear(); }
  void LoadInteractors(std::vector<Object*> ix_vec);
  double const GetCellLength() { return cell_length_; }
  std::vector<pair_interaction> GetPairInteractions();
};

#endif  // _SIMCORE_CELL_LIST_H_
