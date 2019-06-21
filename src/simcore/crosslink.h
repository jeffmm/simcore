#ifndef _SIMCORE_CROSSLINK_H_
#define _SIMCORE_CROSSLINK_H_

//#include "species.h"
#include "anchor.h"
#include "minimum_distance.h"
//#include <mutex>

enum xstate {
  unbound,
  singly,
  doubly
};

class Crosslink : public Object {
  private:
    Interaction ix;
    MinimumDistance * mindist_;
    bool doubly_bound_;
    draw_type draw_;
    bind_state state_;
    double k_on_,
           k_off_,
           k_spring_,
           k_align_,
           f_spring_max_,
           rest_length_,
           rcapture_;
    std::vector<Anchor> anchors_;
    std::vector<Object*> nlist_;
    void CalculateTetherForces();
    void AttemptCrosslink();
    /* TODO, get rid of racy neighborlist additions */
    //std::mutex xlink_mtx_;

  public:
    Crosslink();
    void Init(MinimumDistance * mindist);
    void UnbindAnchor(bool second=false);
    void AttachObjRandom(Object * obj);
    void UpdateCrosslink();
    Anchor * GetBoundPtr();
    void Draw(std::vector<graph_struct*> * graph_array);
    void SetDoubly();
    void SetSingly();
    void SetUnbound();
    bool IsDoubly();
    void UpdatePosition();
    void AddNeighbor(Object * neighbor);
    bind_state GetState();
};

#endif
