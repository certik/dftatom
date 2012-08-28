include_directories(${PROJECT_BINARY_DIR}/src)

project(oscillator)

add_executable(osc_nonrel osc_nonrel.f90)
target_link_libraries(osc_nonrel dftatom)

add_executable(osc_rel osc_rel.f90)
target_link_libraries(osc_rel dftatom)

add_test(osc_nonrel ${PROJECT_BINARY_DIR}/osc_nonrel)
add_test(osc_rel ${PROJECT_BINARY_DIR}/osc_rel)
