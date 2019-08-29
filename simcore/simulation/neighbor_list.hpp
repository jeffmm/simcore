#ifndef _SIMCORE_NEIGHBOR_LIST_H_
#define _SIMCORE_NEIGHBOR_LIST_H_

#include <algorithm>
#include <map>
#include <unordered_set>
#include <vector>

typedef std::unordered_set<int> registry;

typedef std::map<int, registry> verlet_list;

class NeighborList {
 private:
  std::mutex mtx_;
  registry reg_;  // list of object ids in our verlet list
  verlet_list table_;
  verlet_list table_persistent_;
  void Register(int id) {
    std::lock_guard<std::mutex> lk(mtx_);
    reg_.insert(id);
    registry neighbors;
    table_[id] = neighbors;
    table_persistent_[id] = neighbors;
  }
  // returns true if not previously registered (id does not exist in table)
  bool NotRegistered(int id) { return !(reg_.count(id) > 0); }
  bool NotNeighbors(int id1, int id2) { return !(table_[id1].count(id2) > 0); }

 public:
  //  Returns true if neighbors, false if not neighbors
  bool AreNeighbors(int id1, int id2) {
    if (NotRegistered(id1) || NotRegistered(id2))
      return false;
    else if (NotNeighbors(id1, id2))
      return false;
    else
      return true;
  }
  void AddNeighbors(int id1, int id2) {
    if (NotRegistered(id1)) Register(id1);
    if (NotRegistered(id2)) Register(id2);
    std::lock_guard<std::mutex> lk(mtx_);
    table_[id1].insert(id2);
    table_[id2].insert(id1);
  }
  void RemoveNeighbors(int id1, int id2) {
    std::lock_guard<std::mutex> lk(mtx_);
    table_[id1].erase(id2);
    table_[id1].erase(id2);
  }
  void CompareLists(int* inits, int* completes) {
    for (auto id = reg_.begin(); id != reg_.end(); ++id) {
      // First check for mesh ids that are no longer crossing
      for (auto it = table_persistent_[*id].begin();
           it != table_persistent_[*id].end();) {
        if (NotNeighbors(*id, *it)) {
          (*completes)++;
          it = table_persistent_[*id].erase(it);
        } else {
          ++it;
        }
      }
      // Now check for new crossings
      int size = table_persistent_[*id].size();
      for (auto it = table_[*id].begin(); it != table_[*id].end(); ++it) {
        table_persistent_[*id].insert(*it);
      }
      (*inits) += (table_persistent_[*id].size() - size);
      table_[*id].clear();
    }
  }
};

#endif
