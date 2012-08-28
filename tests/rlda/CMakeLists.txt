include_directories(${PROJECT_BINARY_DIR}/src)

project(nist)

add_executable(conv_rlda conv_rlda.f90)
target_link_libraries(conv_rlda dftatom)

add_test(conv_rlda ${PROJECT_BINARY_DIR}/conv_rlda)
