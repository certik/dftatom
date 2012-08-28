include_directories(${PROJECT_BINARY_DIR}/src)

project(hydrogen)
add_executable(hydrogen hydrogen.f90)
target_link_libraries(hydrogen dftatom)

add_test(hydrogen ${PROJECT_BINARY_DIR}/hydrogen)
