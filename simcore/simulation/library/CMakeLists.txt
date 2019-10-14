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
set(LIB ${LIB} GLEW::GLEW ${GLFW3_LIBRARIES} ${OPENGL_gl_LIBRARY} glfw)
set(INCLUDES ${INCLUDES} ${GLEW_INCLUDE_DIRS} ${GLFW3_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS})
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

if(UNIT_TESTS)
    add_definitions(-DUNIT_TESTS=TRUE)
    add_definitions(-DTRACE=TRUE)
    add_definitions(-DDEBUG=TRUE)
endif()

include_directories(${INCLUDES})
set(TARGET "simcore")
set(SOURCES anchor.cpp 
            bond.cpp
            br_bead.cpp
            br_bead_species.cpp
            br_rod.cpp
            cell.cpp
            cell_list.cpp
            cpu.cpp
            cpu_time.cpp
            crosslink.cpp
            crosslink_species.cpp
            crosslink_manager.cpp
            filament.cpp
            filament_species.cpp
            grabber.cpp
            graphics.cpp
            interaction_manager.cpp
            linear_algebra.cpp
            logger.cpp
            mesh.cpp
            minimum_distance.cpp
            object.cpp
            params_parser.cpp
            rng.cpp
            simulation.cpp
            simulation_manager.cpp
            site.cpp
            space.cpp
            species_base.cpp
            spherocylinder.cpp
            spherocylinder_species.cpp
            spindle.cpp
            spindle_species.cpp
            struct_analysis.cpp
            writebmp.cpp
)

add_library(${TARGET} STATIC ${SOURCES})
target_link_libraries(${TARGET} ${LIB})
target_compile_options(${TARGET} PRIVATE -Wno-parentheses -Wno-return-type -Wno-deprecated-declarations)
target_include_directories(${TARGET} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/${TARGET}>
)
