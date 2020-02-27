#ifndef _SIMCORE_CELL_LIST_H_
#define _SIMCORE_CELL_LIST_H_

#include "cell.hpp"
#include "logger.hpp"

class CellList {
private:
  Cell ***cell_;
  static double _min_cell_length_;
  static int _n_dim_;
  static int _n_periodic_;
  static int _n_cells_1d_;
  static double _cell_length_;
  void AllocateCells();
  void DeallocateCells();
  void LabelCells();
  void AssignCellNeighbors(bool redundancy=false);
  xyz_coord FindCellCoords(Object &obj);
  void ClearCellNeighbors();

public:
  CellList() {}
  static void Init(int n_dim, int n_periodic, double system_radius);
  static void SetMinCellLength(double l);
  static double GetCellLength();
  void MakePairs(std::vector<Interaction> &pair_list);
  void RenewObjectsCells(std::vector<Object *> &objs);
  void ResetNeighbors();
  void AssignObjectsCells(std::vector<Object *> &objs);
  void PairSingleObject(Object &obj, std::vector<Interaction> &pair_list);
  void ClearCellObjects();
  void Clear();
  void BuildCellList();
};

#endif
