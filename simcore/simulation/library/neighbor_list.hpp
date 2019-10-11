#include "auxiliary.hpp"
#include "object.hpp"
#include <mutex>

/* A data structure that is used to hold a list of particles that are nearby the
 * owner of the list in a thread-safe way. */
class NeighborList {
private:
  std::mutex mtx_;
  std::vector<Object *> nlist_;

public:
  NeighborList() {}
  ~NeighborList() { Clear(); }
  NeighborList(const NeighborList &that) { this->nlist_ = that.nlist_; }
  NeighborList &operator=(NeighborList const &that) {
    this->nlist_ = that.nlist_;
    return *this;
  }
  void AddNeighbor(Object *obj) {
    std::lock_guard<std::mutex> lk(mtx_);
    nlist_.push_back(obj);
  }
  const Object *const *GetNeighborListMem() { return &nlist_[0]; }
  void Clear() { nlist_.clear(); }
  const int NNeighbors() const { return nlist_.size(); }
  Object *GetNeighbor(int i_neighbor) {
    if (i_neighbor >= nlist_.size()) {
      Logger::Error("Invalid index received in class NeighborList");
    }
    return nlist_[i_neighbor];
  }
};
