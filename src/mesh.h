#ifndef _SIMCORE_MESH_H_
#define _SIMCORE_MESH_H_

#include <vector>

class Object {
  private:
    int static next_oid_;
  protected:
    int oid_;
    double position_[3];
    double orientation_[3];
    double rotation_[3];
    double length_, diameter_, area_;
    void SetOID() { oid_ = next_oid_++; }
  public:
    int const GetOID() { return oid_; }
    void SetPosition(double* pos) {
      std::copy(pos,pos+3,position_);
    }
    void SetOrientation(double *u) {
      std::copy(u, u+3, orientation_);
    }
    void SetRotation(double *r) {
      std::copy(r, r+3, rotation_);
    }
    void SetLength(double l) {
      length_ = l;
    }
    void SetDiameter(double d) {
      diameter_ = d;
    }
    void SetArea(double a) {
      area_ = a;
    }
    Object() {SetOID();}
};

int Object::next_oid_ = 0;
class Bond;
class Triangle;
class Site : public Object {
  // Vertices don't know how bonds or triangles work
  protected:
    std::vector<Bond*> bond_;
    std::vector<Triangle*> triangle_;
  public:
    Site() {}
    void AddBond(Bond *b);
    int GetAdjacentBondOID();
};

class Bond : public Object {
  // Bonds know vertices, but don't know how triangles work
  protected:
    std::vector<Site*> site_;
    std::vector<Triangle*> triangle_;
  public:
    Bond() {}
    void AddSite(Site s) { site_.push_back(&s); }
};

void Site::AddBond(Bond *b) {
  bond_.push_back(b);
}

int Site::GetAdjacentBondOID() {
  int oid=0;
  for (std::vector<Bond*>::iterator it=bond_.begin(); it!=bond_.end(); ++it) {
    oid = (*it)->GetOID();
  }
  return oid;
}

class Triangle : public Object {
  // Triangles know how both bonds and vertices work
  protected:
    std::vector<Bond*> bond_;
    std::vector<Site*> site_;
  public:
    Triangle() {}
};

class Mesh : public Object {
  protected:
    bool filled_; // If true, simples are triangles, else, simples are bonds
    std::vector<Site> site_;
    std::vector<Bond> bond_;
    std::vector<Triangle> triangle_;
  public:
    Mesh() {}
    void AddSite(Site s) {
      site_.push_back(s);
    }
    void AddBond(Bond b) {
      bond_.push_back(b);
    }
    void AddTriangle(Triangle t) {
      triangle_.push_back(t);
    }
};


#endif // _SIMCORE_MESH_H_
