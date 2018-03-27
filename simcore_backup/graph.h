#ifndef _SIMCORE_NETWORK_H_
#define _SIMCORE_NETWORK_H_

class Network {
  protected:
    std::vector<Site> sites_;
    std::vector<Bond> bonds_;

  public:
    Network() {}

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

#endif // _SIMCORE_NETWORK_H_
