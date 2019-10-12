set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(yaml-cpp REQUIRED)
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)

set(LIB ${YAML_CPP_LIBRARIES} ${GSL_LIBRARIES}  ${FFTW_LIBRARIES} KMC)

set(INCLUDES ${GSL_INCLUDE_DIRS} ${FFTW_INCLUDE_DIR} ${YAML_CPP_INCLUDE_DIR})

# If we aren't compiling on local machines, use openmp
if (OMP)
    find_package(OpenMP REQUIRED)
    set(LIB ${LIB} OpenMP::OpenMP_CXX)
endif()

if (GRAPH)
  find_package(glfw3 REQUIRED)
  find_package(glew REQUIRED)
  find_package(OpenGL REQUIRED)
  set(LIB ${LIB} GLEW::GLEW ${GLFW3_LIBRARIES} ${OPENGL_gl_LIBRARY} glfw )
set(INCLUDES ${INCLUDES} ${GLEW_INCLUDE_DIRS} ${GLFW3_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS})
else()
  add_definitions(-DNOGRAPH=TRUE)
endif()
include_directories(${INCLUDES})

set(TARGET simcore.exe)
add_executable(${TARGET} simcore_main.cpp)
target_link_libraries(${TARGET} simcore)
