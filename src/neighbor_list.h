#ifndef _SIMCORE_NEIGHBOR_LIST_H_
#define _SIMCORE_NEIGHBOR_LIST_H_

#include <unordered_set>
#include <map>
#include <vector>
#include <algorithm>

typedef std::unordered_set<int> registry;

typedef std::map<int, registry> verlet_list;

class NeighborList {
  private:
    registry reg_; // list of object ids in our verlet list
    verlet_list table_;
    void Register(int id) {
      reg_.insert(id);
      registry neighbors;
      table_[id] = neighbors;
    }
    // returns true if not previously registered (id does not exist in table)
    bool NotRegistered(int id) {
      if (reg_.find(id) == reg_.end()) return true;
      return false;
    }

  public:
    //  Returns true if neighbors, false if not neighbors
    bool AreNeighbors(int id1, int id2) {
      if (NotRegistered(id1) || NotRegistered(id2)) return false;
      else if (table_[id1].find(id2) == table_[id1].end()) return false;
      else return true;
    }
    void AddNeighbors(int id1, int id2) {
      if (NotRegistered(id1)) Register(id1);
      if (NotRegistered(id2)) Register(id2);
      table_[id1].insert(id2);
      table_[id2].insert(id1);
    }
    void RemoveNeighbors(int id1, int id2) {
      table_[id1].erase(id2);
      table_[id2].erase(id1);
    }
};

#endif
