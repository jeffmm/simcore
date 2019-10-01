#include "cell_list.hpp"

void CellList::Init(int n_cells_1d, double cell_length, int n_dim,
                    int n_periodic) {
  Logger::Trace("Initializing cell list");
  n_dim_ = n_dim;
  n_periodic_ = n_periodic;
  cell_length_ = cell_length;
  n_cells_1d_ = n_cells_1d;
  AllocateCells();
  LabelCells();
  ClearCellObjects();
  /* Build with redundant neighbor pairs for fast overlap checking of new
     objects added to cell list */
  AssignCellNeighbors(true);
}

void CellList::AllocateCells() {
  Logger::Debug("Allocating cell list");
  Logger::Trace("cell_length: %2.2f", cell_length_);
  Logger::Trace("n_cells_1d: %d", n_cells_1d_);

  int third_dim = (n_dim_ == 3 ? n_cells_1d_ : 1);
  cell_ = new Cell **[n_cells_1d_];
  for (int i = 0; i < n_cells_1d_; ++i) {
    cell_[i] = new Cell *[n_cells_1d_];
    for (int j = 0; j < n_cells_1d_; ++j) {
      cell_[i][j] = new Cell[third_dim];
    }
  }
}

void CellList::Clear() {
  ClearCellObjects();
  ClearCellNeighbors();
  DeallocateCells();
}

void CellList::ResetNeighbors() {
  /* Rebuild cell neighbors without redundant neighbor pairs */
  ClearCellNeighbors();
  AssignCellNeighbors(false);
}

void CellList::DeallocateCells() {
  Logger::Debug("Deallocating cell list");
  for (int i = 0; i < n_cells_1d_; ++i) {
    for (int j = 0; j < n_cells_1d_; ++j) {
      delete[] cell_[i][j];
    }
    delete[] cell_[i];
  }
  delete[] cell_;
}

void CellList::LabelCells() {
  int third_dim = (n_dim_ == 3 ? n_cells_1d_ : 1);
  for (int i = 0; i < n_cells_1d_; ++i) {
    for (int j = 0; j < n_cells_1d_; ++j) {
      for (int k = 0; k < third_dim; ++k) {
        cell_[i][j][k].AssignIndex(i, j, k);
      }
    }
  }
}

void CellList::MakePairs(std::vector<Interaction> &pair_list) {
  Logger::Debug("Constructing object interaction pairs");
  int third_dim = (n_dim_ == 3 ? n_cells_1d_ : 1);
  for (int i = 0; i < n_cells_1d_; ++i) {
    for (int j = 0; j < n_cells_1d_; ++j) {
      for (int k = 0; k < third_dim; ++k) {
        cell_[i][j][k].MakePairs(pair_list);
      }
    }
  }
}

xyz_coord CellList::FindCellCoords(Object &obj) {
  const double *const spos = obj.GetScaledPosition();
  double x = spos[0] + 0.5;
  int xcell = (int)floor(n_cells_1d_ * x);
  if (xcell == n_cells_1d_)
    xcell -= 1;
  double y = spos[1] + 0.5;
  int ycell = (int)floor(n_cells_1d_ * y);
  if (ycell == n_cells_1d_)
    ycell -= 1;
  int zcell = 0;
  if (n_dim_ == 3) {
    double z = spos[2] + 0.5;
    zcell = (int)floor(n_cells_1d_ * z);
    if (zcell == n_cells_1d_)
      zcell -= 1;
  }
  return std::make_tuple(xcell, ycell, zcell);
}

void CellList::ClearCellObjects() {
  Logger::Trace("Clearing cell list objects");
  int third_dim = (n_dim_ == 3 ? n_cells_1d_ : 1);
  for (int i = 0; i < n_cells_1d_; ++i) {
    for (int j = 0; j < n_cells_1d_; ++j) {
      for (int k = 0; k < third_dim; ++k) {
        cell_[i][j][k].ClearObjs();
      }
    }
  }
}

void CellList::ClearCellNeighbors() {
  Logger::Trace("Clearing cell list neighbors");
  int third_dim = (n_dim_ == 3 ? n_cells_1d_ : 1);
  for (int i = 0; i < n_cells_1d_; ++i) {
    for (int j = 0; j < n_cells_1d_; ++j) {
      for (int k = 0; k < third_dim; ++k) {
        cell_[i][j][k].ClearNeighbors();
      }
    }
  }
}

void CellList::RenewObjectsCells(std::vector<Object *> &objs) {
  ClearCellObjects();
  AssignObjectsCells(objs);
}

void CellList::AssignObjectsCells(std::vector<Object *> &objs) {
  Logger::Debug("Assigning objects to cells");
  for (auto obj = objs.begin(); obj != objs.end(); ++obj) {
    int x, y, z;
    std::tie(x, y, z) = FindCellCoords(**obj);
#ifdef TRACE
    Logger::Trace("Object %d assigned to %s", (*obj)->GetOID(),
                  cell_[x][y][z].Report().c_str());
#endif
    cell_[x][y][z].AddObj(**obj);
  }
}

/* In order to have a cell list with unique cell neighbors, we want:
   For a given cell at position x_i, y_i, z_i, we want all 9 cells that are
   adjacent to us at (*, *, z_{i+1}), the three y cells adjacent to us at
   (*, y_{i+1}, z_i), and one cell adjacent to us along x at (x_{i+1}, y_i, z_i)
 */
void CellList::AssignCellNeighbors(bool redundancy) {
  /* If we build with redundancy, all cells will store neighbors of all
     adjacent cells, so total neighbor pairs stored will be doubled. This may be
     useful for quickly determining the potential interactions from a single
     object (ie useful for quick overlap checking of new objects) */
  if (redundancy) {
    Logger::Debug(
        "Assigning cell list neighbors with redundant neighbor pairs");
  } else {
    Logger::Debug("Assigning cell list neighbors");
  }
  int third_dim = (n_dim_ == 3 ? n_cells_1d_ : 1);
  // Loop through all cells in cell list
  for (int z = 0; z < third_dim; ++z) {
    for (int y = 0; y < n_cells_1d_; ++y) {
      for (int x = 0; x < n_cells_1d_; ++x) {
        Cell &c = cell_[x][y][z];
        int z_begin = (redundancy && n_dim_ == 3 ? z - 1 : z);
        int z_end = (n_dim_ == 3 ? z + 2 : z + 1);
        /* Add all adjacent cells "above" this cell along z axis */
        for (int zp = z_begin; zp < z_end; ++zp) {
          int nz = zp;
          if ((nz < 0 || nz == n_cells_1d_) && n_periodic_ >= 3) {
            nz = (nz < 0 ? n_cells_1d_ - 1 : 0);
          } else if (nz < 0 || nz == n_cells_1d_) {
            continue;
          }
          int y_begin = (redundancy || zp > z ? y - 1 : y);
          for (int yp = y_begin; yp < y + 2; ++yp) {
            int ny = yp;
            if ((ny < 0 || ny == n_cells_1d_) && n_periodic_ >= 2) {
              ny = (ny < 0 ? n_cells_1d_ - 1 : 0);
            } else if (ny < 0 || ny == n_cells_1d_) {
              continue;
            }
            int x_begin = (redundancy || yp > y || zp > z ? x - 1 : x + 1);
            for (int xp = x_begin; xp < x + 2; ++xp) {
              int nx = xp;
              if ((nx < 0 || nx == n_cells_1d_) && n_periodic_ >= 1) {
                nx = (nx < 0 ? n_cells_1d_ - 1 : 0);
              } else if (nx < 0 || nx == n_cells_1d_) {
                continue;
              }
              /* Note that cells never add themselves */
              if (nx == x && ny == y && nz == z) {
                continue;
              }
              c.AddNeighbor(cell_[nx][ny][nz]);
              Logger::Trace("%s has neighbor %s", c.Report().c_str(),
                            cell_[nx][ny][nz].Report().c_str());
            }
          }
        }
      }
    }
  }
}

/* Relies on redundant cell list neighbor pairs to identify all potential
   interactions with this object */
void CellList::PairSingleObject(Object &obj,
                                std::vector<Interaction> &pair_list) {
  int x, y, z;
  std::tie(x, y, z) = FindCellCoords(obj);
  Logger::Trace("Making pairs with single object %d in %s", obj.GetOID(),
                cell_[x][y][z].Report().c_str());
  const std::vector<Cell *> neighbors = cell_[x][y][z].GetCellNeighbors();
  cell_[x][y][z].PairSingleObject(obj, pair_list);
  for (auto cell = neighbors.begin(); cell != neighbors.end(); ++cell) {
    (*cell)->PairSingleObject(obj, pair_list);
  }
}
