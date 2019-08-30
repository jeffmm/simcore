#include "cell_list.hpp"

Cell::Cell() { objs_.reserve(5); }

void Cell::PairObjs(Cell* c, std::vector<ix_pair>* nl) {
  int n_objs_ = objs_.size();
  if (n_objs_ == 0) return;
  if (c == this) {
    for (int i = 0; i < n_objs_; ++i) {
      for (int j = i + 1; j < n_objs_; ++j) {
        nl->push_back(std::make_pair(objs_[i], objs_[j]));
      }
    }
    return;
  }
  int m_objs_ = c->objs_.size();
  if (m_objs_ == 0) return;
  for (int i = 0; i < n_objs_; ++i) {
    for (int j = 0; j < m_objs_; ++j) {
      nl->push_back(std::make_pair(objs_[i], c->objs_[j]));
    }
  }
}

void Cell::AddObj(int oid) {
  std::lock_guard<std::mutex> lk(cell_mtx);
  objs_.push_back(oid);
}

void Cell::PopBack() { objs_.pop_back(); }
