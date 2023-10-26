include_directories(${PROJECT_BINARY_DIR}/src)

project(well)

set(SRC
    special.f90 optimize.f90
    )

if (CMAKE_Fortran_COMPILER_ID STREQUAL "LFortran")
    add_definitions("--implicit-typing")
endif()

add_library(well ${SRC})
add_dependencies(well dftatom)

add_executable(well_nonrel well_nonrel.f90)
target_link_libraries(well_nonrel dftatom well)

add_executable(well_rel well_rel.f90)
target_link_libraries(well_rel dftatom well)

add_test(well_nonrel ${PROJECT_BINARY_DIR}/well_nonrel)
add_test(well_rel ${PROJECT_BINARY_DIR}/well_rel)
