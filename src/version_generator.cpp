#include "simcore/version_generator.hpp"
/* This file should be configured by CMake to generate a cpp file that gets
   compiled to make the simcore_version variable (the most recent git commit
   sha1) available to any file that includes the above header */
#define GIT_SHA1 "d5dd0b1dfa12be1343316b257c47c9c63075d8e2"
const char simcore_version[] = GIT_SHA1;
