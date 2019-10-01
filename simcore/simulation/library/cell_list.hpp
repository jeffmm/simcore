#ifndef _SIMCORE_CELL_LIST_H_
#define _SIMCORE_CELL_LIST_H_

#include "cell.hpp"

class CellList {
private:
  Cell ***cell_;
  int n_dim_;
  int n_periodic_;
  int n_cells_1d_;
  double cell_length_;
  void AllocateCells();
  void DeallocateCells();
  void LabelCells();
  void AssignCellNeighbors(bool redundancy=false);
  xyz_coord FindCellCoords(Object &obj);
  void ClearCellNeighbors();

public:
  CellList() {}
  void Init(int n_cells_1d, double cell_length, int n_dim, int n_periodic);
  void MakePairs(std::vector<Interaction> &pair_list);
  void RenewObjectsCells(std::vector<Object *> &objs);
  void ResetNeighbors();
  void AssignObjectsCells(std::vector<Object *> &objs);
  void PairSingleObject(Object &obj, std::vector<Interaction> &pair_list);
  void ClearCellObjects();
  void Clear();
};

#endif
