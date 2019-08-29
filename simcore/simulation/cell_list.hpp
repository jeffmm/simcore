#include <mutex>
#include <tuple>
#include <vector>

typedef std::pair<int, int> ix_pair;

class Cell {
 private:
  std::mutex cell_mtx;
  void PairSelf(std::vector<ix_pair>& nl);

 public:
  std::vector<int> objs_;
  Cell();
  void AddObj(int oid);
  void PopBack();
  void PairObjs(Cell* c, std::vector<ix_pair>* nl);
};
