
#ifndef _CYTOSCORE_SIMPLE_OBJECT_H_
#define _CYTOSCORE_SIMPLE_OBJECT_H_

template <typename T>


class Simple {
  protected:
    double pos_[2];
    double force_[2];
  public:
    UpdatePosition();
}

class SpringParticle : public Simple {
  private:
    using Simple::pos_;
    using Simple::force_;
    double mass_;
    double eq_pos_[2];
  public:
    using Simple::UpdatePosition();
}

template <typename T, template<typename> class V>
class Composite {
  protected:
    std::vector<T> parts_;
    V<T> integrator_;
}



#endif
