include_directories(${PROJECT_BINARY_DIR}/src)

project(double_min)

add_executable(rmin_conv rmin_conv.f90)
target_link_libraries(rmin_conv dftatom)
add_test(rmin_conv ${PROJECT_BINARY_DIR}/rmin_conv)

add_executable(N_conv N_conv.f90)
target_link_libraries(N_conv dftatom)
add_test(N_conv ${PROJECT_BINARY_DIR}/N_conv)

add_executable(fd fd.f90)
target_link_libraries(fd dftatom)
add_test(fd ${PROJECT_BINARY_DIR}/fd)
