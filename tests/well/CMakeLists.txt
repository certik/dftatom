include_directories(${PROJECT_BINARY_DIR}/src)

project(well)

add_executable(well_nonrel well_nonrel.f90 special.f90 optimize.f90)
target_link_libraries(well_nonrel dftatom)

add_executable(well_rel well_rel.f90 special.f90 optimize.f90)
target_link_libraries(well_rel dftatom)

add_test(well_nonrel ${PROJECT_BINARY_DIR}/well_nonrel)
add_test(well_rel ${PROJECT_BINARY_DIR}/well_rel)
