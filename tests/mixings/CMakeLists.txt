include_directories(${PROJECT_BINARY_DIR}/src)

project(mixings)

add_executable(mixings mixings.f90)
target_link_libraries(mixings dftatom)

add_test(mixings ${PROJECT_BINARY_DIR}/mixings)
