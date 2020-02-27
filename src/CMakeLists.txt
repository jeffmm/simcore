find_package(yaml-cpp REQUIRED)
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)

set(LIB ${YAML_CPP_LIBRARIES} ${GSL_LIBRARIES}  ${FFTW_LIBRARIES} KMC)
set(INCLUDES ${GSL_INCLUDE_DIRS} ${FFTW_INCLUDE_DIR} ${YAML_CPP_INCLUDE_DIR})

if (OMP)
    find_package(OpenMP REQUIRED)
    set(LIB ${LIB} OpenMP::OpenMP_CXX)
endif()

if (GRAPH)
  find_package(glfw3 REQUIRED)# HINTS /usr/include/GLFW)
  find_package(OpenGL REQUIRED)
  if (WINDOWS_MODE)
    message("Using Windows Graphics")
    set(GLEW_LIBRARY -lGL -lGLU -lGLEW -ldl -lX11 -lpthread -lXrandr -lXxf86vm -lXi -lm -lXcursor -lXinerama)
    set(GLEW_INCLUDE_PATH /usr/include/)
    add_definitions(-DWINGRAPH=TRUE)
  else()
    find_package(glew REQUIRED)
  endif()
  set(LIB ${LIB} ${GLEW_LIBRARY} ${GLFW3_LIBRARIES} ${OPENGL_gl_LIBRARY} glfw)
  set(INCLUDES ${INCLUDES} ${GLEW_INCLUDE_PATH} ${GLFW3_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS})
else()
  add_definitions(-DNOGRAPH=TRUE)
endif()

if(DEBUG)
    add_definitions(-DDEBUG=TRUE)
endif()

if(TRACE)
    add_definitions(-DTRACE=TRUE)
    add_definitions(-DDEBUG=TRUE)
endif()

if (TESTS)
    add_definitions(-DTESTS=TRUE)
endif()

# Generate version definition from git sha1. This
# needs to be configured before the cpp files are
# GLOB'd into the SOURCES variable or explicitly
# appended to the SOURCES list
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/version_generator.cpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/version_generator.cpp" @ONLY)

# Find all headers in the project include folder
file(GLOB HEADERS CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/simcore/*.hpp")
file(GLOB SOURCES "*.cpp")

# Make a static library (tradeoff of memory for better performance)
add_library(simcore STATIC ${SOURCES} ${HEADERS})

# Add an alias so the target can be used while building
add_library(simcore::simcore ALIAS simcore)

target_link_libraries(simcore ${LIB})

# Set the include directory for builds (tests, examples, etc) and for clients
target_include_directories(simcore PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Associate the include files with the target, so they will be installed
set_target_properties(simcore PROPERTIES PUBLIC_HEADER "${HEADERS}")

# Set install targets and destinations
install(TARGETS simcore
    EXPORT simcoreTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include/simcore)

# Export additional install files (for clients using CMake)
install(EXPORT simcoreTargets
    NAMESPACE simcore::
    FILE simcoreConfig.cmake
    DESTINATION lib/cmake/simcore)

# Add the executable
add_subdirectory("executable")

# Add the configurator utility
add_subdirectory("configurator")
