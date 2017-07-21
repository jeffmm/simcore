#include "mesh.h"
#include <iostream>

int main() {
  Mesh m;
  Object obj;
  Site s;
  Bond b;
  std::cout << obj.GetOID() << " " << s.GetOID() << "\n";

  b.AddSite(s);
  std::cout << b.GetOID() << "\n";
  s.AddBond(&b);
  m.AddSite(s);
  std::cout << s.GetAdjacentBondOID() << "\n";
  m.AddBond(b);
  return 0;

}
