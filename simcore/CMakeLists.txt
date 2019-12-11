if(CMAKE_BUILD_TYPE MATCHES "^[Dd]ebug")
  message("Building simcore in Debug mode")
else()
  # otherwise set build to release mode by default
  message("Building simcore in Release mode")
endif()

add_subdirectory(configurator)
add_subdirectory(simulation)
