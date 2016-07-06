#include "cell_list.h"

int const Cell::GetNInteractions() {
  int count = 0;
  for (auto it=simples_.begin(); it!=simples_.end();) {
    count += Count()-1;
    for (auto it=neighbors_.begin(); it!=neighbors_.end(); ++it)
      count += (*it)->Count();
    simples_.erase(it);
  }
  return count;
}
std::vector<cell_interaction> Cell::PairInteractions() {
  std::vector<cell_interaction> interactions;
  for (auto it=simples_.begin(); it!=simples_.end();) {
    for (auto jt=it+1; jt!=simples_.end(); ++jt) {
      cell_interaction cell_int(*it,*jt);
      interactions.push_back(cell_int);
    }
    for (auto nt=neighbors_.begin(); nt!=neighbors_.end(); ++nt) {
      for (auto jt=(*nt)->Begin(); jt!=(*nt)->End(); ++jt) {
        cell_interaction cell_int(*it, *jt);
        interactions.push_back(cell_int);
      }
    }
    simples_.erase(it);
  }
  return interactions;
}

void CellList::SetNeighbors(int i, int j, int k) {
  cell_index ind = {i,j,k};
  for (int l=i-1; l<i+2; ++l) {
    int temp_l=l;
    /* Correct for periodicity, or skip cell if out of range */
    if (l<0 && n_periodic_ > 0)
      temp_l = n_cells_1d_-1;
    else if (l>n_cells_1d_-1 && n_periodic_ > 0)
      temp_l = 0;
    else if (l<0 || l>n_cells_1d_-1) continue;
    for (int m=j-1; m<j+2; ++m) {
      int temp_m=m;
      /* Correct for periodicity, or skip cell if out of range*/
      if (m<0 && n_periodic_ > 1)
        temp_m = n_cells_1d_-1;
      else if (m>n_cells_1d_-1 && n_periodic_ > 1)
        temp_m = 0;
      else if (m<0 || m>n_cells_1d_-1) continue;
      if (n_dim_==2) {
        if (l==i && m==j) continue;
        cell_index neighbor_ind = {temp_l,temp_m,0};
        Cell * c = cells_[neighbor_ind].GetCellPtr();
        cells_[ind].AddNeighbor(c);
      }
      else if (n_dim_==3) {
        for (int n=k-1; n<k+2; ++n) {
          int temp_n = n;
          /* Correct for periodicity, or skip if out of range */
          if (n<0 && n_periodic_ > 2)
            temp_n = n_cells_1d_-1;
          else if (n>n_cells_1d_-1 && n_periodic_ > 1)
            temp_n = 0;
          else if (n<0 || n>n_cells_1d_-1) continue;
          /* Don't add self to neighbor list */
          if (l==i && m==j && n==k) continue;
          /* Add pointer to neighbor to cell neighbor list */
          cell_index neighbor_ind = {temp_l,temp_m,temp_n};
          Cell * c = (cells_[neighbor_ind]).GetCellPtr();
          cells_[ind].AddNeighbor(c);
        }
      }
    }
  }
}

void CellList::Init(int n_dim, int n_periodic, double cell_length, double system_radius) {
  n_dim_ = n_dim;
  n_periodic_ = n_periodic;
  n_cells_1d_ = (int) ceil(2.0*system_radius/cell_length);
  if (n_cells_1d_ < 3)
    n_cells_1d_ = 3;
  cell_length_ = 2.0*system_radius/n_cells_1d_;
  scaled_cell_length_ = 1.0/n_cells_1d_;
  ClearCells();
  InitCells();
}

void CellList::InitCells() {
  for (int i=0; i<n_cells_1d_; ++i) {
    for (int j=0; j<n_cells_1d_; ++j) {
      if (n_dim_==2) {
        cell_index ind = {i,j,0};
        Cell c(i,j,0);
        cells_[ind] = c;
      }
      else if (n_dim_==3) {
        for (int k=0; k<n_cells_1d_; ++k) {
          cell_index ind = {i,j,k};
          Cell c(i,j,k);
          cells_[ind] = c;
        }
      }
    }
  }
  /* Set the nearest neighbors only after all
     cells are initialized */
  for (int i=0; i<n_cells_1d_; ++i) {
    for (int j=0; j<n_cells_1d_; ++j) {
      if (n_dim_==2) {
        SetNeighbors(i,j,0);
      }
      else if (n_dim_==3) {
        for (int k=0; k<n_cells_1d_; ++k)
          SetNeighbors(i,j,k);
      }
    }
  }
}

void CellList::LoadSimples(std::vector<Simple*> sim_vec) {
  for (auto sim = sim_vec.begin(); sim!= sim_vec.end(); ++sim) {
    double const * const s = (*sim)->GetScaledPosition();
    cell_index ind = {0,0,0};
    for (int i=0; i<n_dim_; ++i)
      ind[i]=(int) floor((s[i]+0.5)/scaled_cell_length_);
    cells_[ind].AddSimple(*sim);
  }
}

std::vector<cell_interaction> CellList::GetInteractions() {
  std::vector<cell_interaction> interactions;
  for (cell_map_it it=cells_.begin(); it!= cells_.end(); ++it) {
    std::vector<cell_interaction> j = it->second.PairInteractions();
    interactions.insert(interactions.end(), j.begin(), j.end());
  }
  return interactions;
}

