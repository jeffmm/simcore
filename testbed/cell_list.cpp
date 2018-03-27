#include "cell_list.h"

Cell::Cell() {
  objs.reserve(5);
}

void Cell::PairObjs(Cell * c, std::vector<ix_pair> * nl) {
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

void Cell::AddObj(int oid) {
  std::lock_guard<std::mutex> lk(cell_mtx);
  objs.push_back(oid);
}
