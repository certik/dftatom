include_directories(${PROJECT_BINARY_DIR}/src)

project(Z92n7all)

add_executable(Z92_nonrel Z92_nonrel.f90)
target_link_libraries(Z92_nonrel dftatom)

add_executable(Z92_nonrel_100 Z92_nonrel_100.f90)
target_link_libraries(Z92_nonrel_100 dftatom)

add_executable(Z92_rel Z92_rel.f90)
target_link_libraries(Z92_rel dftatom)

add_executable(Z92_rel_100 Z92_rel_100.f90)
target_link_libraries(Z92_rel_100 dftatom)

add_executable(Z92_tf Z92_tf.f90)
target_link_libraries(Z92_tf dftatom)

add_test(Z92_nonrel ${PROJECT_BINARY_DIR}/Z92_nonrel)
add_test(Z92_tf ${PROJECT_BINARY_DIR}/Z92_tf)
add_test(Z92_rel ${PROJECT_BINARY_DIR}/Z92_rel)
add_test(Z92_nonrel_100 ${PROJECT_BINARY_DIR}/Z92_nonrel_100)
add_test(Z92_rel_100 ${PROJECT_BINARY_DIR}/Z92_rel_100)
