#include "object.h"
#include "omp.h"
#include <iostream>
#include <math.h>
#include <tuple>
#include <mutex>
typedef std::pair<int,int> ix_pair;

class Cell {
  private:
    std::mutex cell_mtx;
    void PairSelf(std::vector<ix_pair> & nl) {
    }

  public:
    std::vector<int> objs;
    Cell() {
      objs.reserve(5);
    }
    void AddObj(int oid) {
      std::lock_guard<std::mutex> lk(cell_mtx);
      objs.push_back(oid);
    }
    void PairObjs(Cell * c, std::vector<ix_pair> * nl) {
      int n_objs = objs.size();
      if (n_objs == 0) return;
      if (c == this) {
        for (int i=0; i<n_objs; ++i) {
          for (int j=i+1; j<n_objs; ++j) {
            nl->push_back(std::make_pair(objs[i],objs[j]));
          }
        }
        return;
      }
      int m_objs = c->objs.size();
      if (m_objs == 0) return;
      for (int i=0; i<n_objs; ++i) {
        for (int j=0; j<m_objs; ++j) {
          nl->push_back(std::make_pair(objs[i],c->objs[j]));
        }
      }
    }
};

class Simulation {
  private:
    std::mutex sim_mtx;
    int count = 0;
    double cell_length_1d;
    int n_cells_1d;
    std::vector<Object> objs;
    std::vector<ix_pair> nlist;
    Cell *** clist;
    void InitRNG();
    void InitPositions();
    void CreatePairs();
    void CreatePairsCellList();
    void AssignCells();
    void PrintPairs();
    void InteractPairs();
    void Interact(int i, int j, int * count);
    void AllocateCellList();
    void DeallocateCellList();
  public:
    void Run();
};
