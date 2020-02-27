# Add the simcore executable
add_executable(simcore.exe "main.cpp")
target_link_libraries(simcore.exe simcore)
install(TARGETS simcore.exe DESTINATION bin)
