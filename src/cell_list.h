#ifndef _SIMCORE_CELL_LIST_H_
#define _SIMCORE_CELL_LIST_H_

#include "auxiliary.h"
#include "object.h"
#include "interaction.h"
//#include <unordered_map>

typedef std::array<int,3> cell_index;
typedef std::pair<Simple*,Simple*> simple_pair;
typedef std::pair< simple_pair , Interaction > pair_interaction;

class Cell {
  private:
    cell_index index_;
    std::vector<Simple*> simples_;
    std::vector<Cell*> neighbors_;
    std::vector< simple_pair > interactions_;
  public:
    Cell() {}
    Cell(int i, int j, int k) {index_ = {i,j,k};}
    void AddSimple(Simple * sim) {simples_.push_back(sim);}
    void AddNeighbor(Cell * neighbor) {neighbors_.push_back(neighbor);}
    std::vector<Simple*>::iterator Begin() {return simples_.begin();}
    std::vector<Simple*>::iterator End() {return simples_.end();}
    int const Count() {return simples_.size();}
    Cell * GetCellPtr() {return this;}
    int const GetNInteractions();
    std::vector<simple_pair> PairInteractions();
};

//typedef std::unordered_map<cell_index, Cell> cell_map;
//typedef std::unordered_map<cell_index, Cell>::iterator cell_map_it;
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
    void Init(int n_dim, int n_periodic, double cell_length, double system_radius);
    void ClearCells() {cells_.clear();}
    void LoadSimples(std::vector<Simple*> sim_vec);
    double const GetCellLength() { return cell_length_; }
    std::vector<pair_interaction> GetPairInteractions();

};

#endif // _SIMCORE_CELL_LIST_H_
