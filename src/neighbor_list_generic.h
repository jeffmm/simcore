// Generic information for neighbor lists
#ifndef _SIMCORE_NEIGHBORLIST_GENERIC_H_
#define _SIMCORE_NEIGHBORLIST_GENERIC_H_

struct _neighbor {
    int idx_;
    double value_; // exists to piggyback calculation onto
};
typedef struct _neighbor neighbor_t;

typedef std::vector<neighbor_t> nl_list;

#endif
