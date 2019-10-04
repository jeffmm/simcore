set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
find_package(yaml-cpp REQUIRED)
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)

set(LIB ${YAML_CPP_LIBRARIES} ${GSL_LIBRARIES}  ${FFTW_LIBRARIES} KMC)

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
endif()
include_directories(${INCLUDES})
set(TARGET "simcore")
set(SOURCES anchor.cpp 
            bead_spring.cpp
            bond.cpp
            br_bead.cpp
            cell.cpp
            cell_list.cpp
            #centrosome.cpp
            cpu.cpp
            cpu_time.cpp
            crosslink.cpp
            crosslink_species.cpp
            crosslink_manager.cpp
            filament.cpp
            filament_species.cpp
            generate_random_unit_vector.cpp
            grabber.cpp
            graphics.cpp
            interaction_engine.cpp
            linear_algebra.cpp
            logger.cpp
            mesh.cpp
            minimum_distance.cpp
            object.cpp
            output_manager.cpp
            params_parser.cpp
            rng.cpp
            simulation.cpp
            simulation_manager.cpp
            site.cpp
            space.cpp
            #species.cpp
            species_base.cpp
            spherocylinder.cpp
            #spindle.cpp
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
