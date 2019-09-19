#include "rng.hpp"

long RNG::_seed_ = 7777777;
std::mutex RNG::_rng_mtx_;
