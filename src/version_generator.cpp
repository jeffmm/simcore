#include "simcore/version_generator.hpp"
/* This file should be configured by CMake to generate a cpp file that gets
   compiled to make the simcore_version variable (the most recent git commit
   sha1) available to any file that includes the above header */
#define GIT_SHA1 "56e69ecb7bcba6dea593ee4e2b8ce5445fb28c43"
const char simcore_version[] = GIT_SHA1;
