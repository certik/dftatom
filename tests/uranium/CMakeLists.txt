include_directories(${PROJECT_BINARY_DIR}/src)

project(uranium)

add_executable(uranium_nonrel uranium_nonrel.f90)
target_link_libraries(uranium_nonrel dftatom)

add_executable(uranium_rel uranium_rel.f90)
target_link_libraries(uranium_rel dftatom)

add_executable(uranium_tf uranium_tf.f90)
target_link_libraries(uranium_tf dftatom)

add_test(uranium_nonrel ${PROJECT_BINARY_DIR}/uranium_nonrel)
add_test(uranium_tf ${PROJECT_BINARY_DIR}/uranium_tf)
add_test(uranium_rel ${PROJECT_BINARY_DIR}/uranium_rel)
