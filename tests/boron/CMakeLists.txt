include_directories(${PROJECT_BINARY_DIR}/src)

project(boron)

add_executable(boron_nonrel boron_nonrel.f90)
target_link_libraries(boron_nonrel dftatom)

add_executable(boron_rel boron_rel.f90)
target_link_libraries(boron_rel dftatom)

add_test(boron_nonrel ${PROJECT_BINARY_DIR}/boron_nonrel)
add_test(boron_rel ${PROJECT_BINARY_DIR}/boron_rel)
