set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
find_package(yaml-cpp REQUIRED)
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)

include_directories(catch2 ${GSL_INCLUDE_DIRS} ${FFTW_INCLUDE_DIR} ${YAML_CPP_INCLUDE_DIR})
set(LIBS ${YAML_CPP_LIBRARIES} ${GSL_LIBRARIES} ${FFTW_LIBRARIES} kmcx)
add_definitions(-DNOGRAPH=TRUE)
add_executable(test_simcore       test_simcore.cpp
                                  ../src/simcore/anchor.cpp 
                                  ../src/simcore/bead_spring.cpp
                                  ../src/simcore/bond.cpp
                                  ../src/simcore/br_bead.cpp
                                  ../src/simcore/cell_list.cpp
                                  ../src/simcore/cpu.cpp
                                  ../src/simcore/cpu_time.cpp
                                  ../src/simcore/crosslink.cpp
                                  ../src/simcore/crosslink_manager.cpp
                                  ../src/simcore/error_exit.cpp
                                  ../src/simcore/filament.cpp
                                  ../src/simcore/filament_species.cpp
                                  ../src/simcore/generate_random_unit_vector.cpp
                                  ../src/simcore/grabber.cpp
                                  ../src/simcore/graphics.cpp
                                  ../src/simcore/interaction_engine.cpp
                                  ../src/simcore/linear_algebra.cpp
                                  ../src/simcore/mesh.cpp
                                  ../src/simcore/minimum_distance.cpp
                                  ../src/simcore/object.cpp
                                  ../src/simcore/output_manager.cpp
                                  ../src/simcore/particle_tracker.cpp
                                  ../src/simcore/passive_filament.cpp
                                  ../src/simcore/rigid_mesh.cpp
                                  ../src/simcore/rng.cpp
                                  #../src/simcore/simcore_main.cpp
                                  ../src/simcore/simulation.cpp
                                  ../src/simcore/simulation_manager.cpp
                                  ../src/simcore/site.cpp
                                  ../src/simcore/space.cpp
                                  ../src/simcore/species.cpp
                                  ../src/simcore/spherocylinder.cpp
                                  ../src/simcore/struct_analysis.cpp
                                  ../src/simcore/writebmp.cpp
)
target_link_libraries(test_simcore ${LIBS})
target_compile_options(test_simcore PRIVATE -Wno-parentheses -Wno-return-type -Wno-deprecated-declarations)
add_test(simcore test_simcore)
