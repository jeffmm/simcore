#include <simcore/cell.hpp>

Cell::Cell() {}

void Cell::AddObj(Object &obj) {
  std::lock_guard<std::mutex> lk(cell_mtx_);
  cell_objs_.push_back(&obj);
}

void Cell::PopBack() { cell_objs_.pop_back(); }

const int Cell::NObjs() const { return cell_objs_.size(); }

void Cell::MakePairs(std::vector<Interaction> &pair_list) const {
  if (NObjs() == 0)
    return;
  MakePairsSelf(pair_list);
  for (auto cell = cell_neighbors_.begin(); cell != cell_neighbors_.end();
       ++cell) {
    MakePairsCell(**cell, pair_list);
  }
}

const std::vector<Cell *> &Cell::GetCellNeighbors() const {
  return cell_neighbors_;
}

const std::vector<Object *> &Cell::GetCellObjects() const { return cell_objs_; }

void Cell::PairSingleObject(Object &obj,
                            std::vector<Interaction> &pair_list) const {
  Logger::Trace("Checking single object pairs with %s", Report().c_str());
  for (int i = 0; i < cell_objs_.size(); ++i) {
    Interaction ix(&obj, cell_objs_[i]);
    pair_list.push_back(ix);
#ifdef TRACE
    Logger::Trace("Single object interaction pair: %d -> %d", obj.GetOID(),
                  cell_objs_[i]->GetOID());
#endif
  }
}

void Cell::MakePairsCell(Cell &cell, std::vector<Interaction> &pair_list) const {
  if (cell.NObjs() == 0)
    return;
  const std::vector<Object *> those_objs = cell.GetCellObjects();
  Logger::Trace("%s adjacent to %s:", Report().c_str(), cell.Report().c_str());
  for (int i = 0; i < NObjs(); ++i) {
    for (int j = 0; j < cell.NObjs(); ++j) {
      Interaction ix(cell_objs_[i], those_objs[j]);
      pair_list.push_back(ix);
#ifdef TRACE
      Logger::Trace("Interaction pair: %d -> %d", cell_objs_[i]->GetOID(),
                    those_objs[j]->GetOID());
#endif
    }
  }
}

void Cell::MakePairsSelf(std::vector<Interaction> &pair_list) const {
  if (NObjs() <= 1) {
    return;
  }
  for (int i = 0; i < NObjs() - 1; ++i) {
    for (int j = i + 1; j < NObjs(); ++j) {
      Interaction ix(cell_objs_[i], cell_objs_[j]);
      pair_list.push_back(ix);
    }
  }
}

void Cell::AssignIndex(const int x, const int y, const int z) {
  cell_index_ = std::make_tuple(x, y, z);
}

std::string Cell::Report() const {
  int x, y, z;
  std::tie(x, y, z) = cell_index_;
  std::string id;
  id = "Cell<" + std::to_string(x) + " " + std::to_string(y) + " " +
       std::to_string(z) + ">";
  return id;
}

void Cell::AddNeighbor(Cell &c) { cell_neighbors_.push_back(&c); }

void Cell::ClearObjs() { cell_objs_.clear(); }
void Cell::ClearNeighbors() { cell_neighbors_.clear(); }
