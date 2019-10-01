#ifndef _SIMCORE_CELL_H_
#define _SIMCORE_CELL_H_

#include "object.hpp"

// ix_pair already defined in interaction.hpp
//typedef std::pair<Object *, Object *> ix_pair;
typedef std::tuple<int, int, int> xyz_coord;

class Cell {
private:
  xyz_coord cell_index_;
  std::mutex cell_mtx_;
  std::vector<Object *> cell_objs_;
  std::vector<Cell *> cell_neighbors_;
  void MakePairsCell(Cell &cell, std::vector<Interaction> &pair_list) const;
  void MakePairsSelf(std::vector<Interaction> &pair_list) const;
  const int NObjs() const;
  const std::vector<Object *> &GetCellObjects() const;

public:
  Cell();
  void AddObj(Object &obj);
  void PopBack();
  void AssignIndex(const int x, const int y, const int z);
  void MakePairs(std::vector<Interaction> &pair_list) const;
  void AddNeighbor(Cell &c);
  std::string Report() const;
  const std::vector<Cell *> &GetCellNeighbors() const;
  void PairSingleObject(Object &obj, std::vector<Interaction> &pair_list) const;
  void ClearObjs();
  void ClearNeighbors();
};

#endif
