#ifndef _SIMCORE_TRACKING_SCHEME_ALLPAIRS_H_
#define _SIMCORE_TRACKING_SCHEME_ALLPAIRS_H_

#include "tracking_scheme.h"

class TrackingSchemeAllPairs : public TrackingScheme {

  public:
    
    TrackingSchemeAllPairs() {
      name_ = "AllPairs";
    }
    virtual ~TrackingSchemeAllPairs() {}

    virtual void Print();
};

#endif
