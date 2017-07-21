#ifndef _SIMCORE_MESH_H_
#define _SIMCORE_MESH_H_

class Mesh {
  protected:
    std::vector<Site> sites_;
    std::vector<Bond> bonds_;

  public:
    Mesh() {}

};

class Bond {
  protected:
    Site *s1_;
    Site *s2_;

  public:
    Bond() {}

};

class Site {
  protected:
    std::vector<Bond> bonds_;

  public:
    Site() {}

};

#endif // _SIMCORE_MESH_H_
