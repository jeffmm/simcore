set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
find_package(yaml-cpp REQUIRED)
find_package(GSL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glew REQUIRED)
find_package(OpenGL REQUIRED)
find_package(FFTW REQUIRED)

include_directories(${GSL_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS} ${GLFW3_INCLUDE_DIRS} ${FFTW_INCLUDE_DIR} ${YAML_CPP_INCLUDE_DIR} ${OPENGL_INCLUDE_DIRS})
set(LIBS ${YAML_CPP_LIBRARIES} ${GSL_LIBRARIES} GLEW::GLEW ${GLFW3_LIBRARIES} ${FFTW_LIBRARIES} ${OPENGL_gl_LIBRARY} glfw)

include_directories(catch2)
add_executable(test_simcore       test_simcore.cpp
                                  ../src/simcore/bead_spring.cpp
                                  ../src/simcore/bond.cpp
                                  ../src/simcore/br_bead.cpp
                                  ../src/simcore/cell_list.cpp
                                  ../src/simcore/centrosome.cpp
                                  ../src/simcore/cpu.cpp
                                  ../src/simcore/cpu_time.cpp
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
                                  ../src/simcore/motor.cpp
                                  ../src/simcore/object.cpp
                                  ../src/simcore/output_manager.cpp
                                  ../src/simcore/particle_tracker.cpp
                                  ../src/simcore/passive_filament.cpp
                                  ../src/simcore/rigid_mesh.cpp
                                  ../src/simcore/simulation.cpp
                                  ../src/simcore/simulation_manager.cpp
                                  ../src/simcore/site.cpp
                                  ../src/simcore/space.cpp
                                  ../src/simcore/species.cpp
                                  ../src/simcore/spherocylinder.cpp
                                  ../src/simcore/spindle.cpp
                                  ../src/simcore/struct_analysis.cpp
                                  ../src/simcore/writebmp.cpp
)
target_link_libraries(test_simcore ${LIBS})
add_test(simcore test_simcore)
