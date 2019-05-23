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
add_executable(simcore  bead_spring.cpp
                        bond.cpp
                        br_bead.cpp
                        cell_list.cpp
                        centrosome.cpp
                        cpu.cpp
                        cpu_time.cpp
                        error_exit.cpp
                        filament.cpp
                        filament_species.cpp
                        generate_random_unit_vector.cpp
                        grabber.cpp
                        graphics.cpp
                        interaction_engine.cpp
                        linear_algebra.cpp
                        mesh.cpp
                        minimum_distance.cpp
                        motor.cpp
                        object.cpp
                        output_manager.cpp
                        particle_tracker.cpp
                        passive_filament.cpp
                        rigid_mesh.cpp
                        simcore_main.cpp
                        simulation.cpp
                        simulation_manager.cpp
                        site.cpp
                        space.cpp
                        species.cpp
                        spherocylinder.cpp
                        spindle.cpp
                        struct_analysis.cpp
                        writebmp.cpp
)
target_link_libraries(simcore ${LIBS})
