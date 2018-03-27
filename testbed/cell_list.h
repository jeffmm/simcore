#include <vector>
#include <mutex>
#include <tuple>

typedef std::pair<int,int> ix_pair;

class Cell {
  private:
    std::mutex cell_mtx;
    void PairSelf(std::vector<ix_pair> & nl);

  public:
    std::vector<int> objs;
    Cell();
    void AddObj(int oid);
    void PairObjs(Cell * c, std::vector<ix_pair> * nl);
};
