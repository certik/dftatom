include_directories(${PROJECT_BINARY_DIR}/src)

project(pseudopotential)

# The erf() function is part of the F2008 standard
add_definitions("-std=f2008")

add_executable(vloc_nonrel vloc_nonrel.f90)
target_link_libraries(vloc_nonrel dftatom)
add_test(vloc_nonrel ${PROJECT_BINARY_DIR}/vloc_nonrel)

add_executable(vloc_nonrel_dft vloc_nonrel_dft.f90)
target_link_libraries(vloc_nonrel_dft dftatom)
add_test(vloc_nonrel_dft ${PROJECT_BINARY_DIR}/vloc_nonrel_dft)

add_executable(vloc_nonrel_dft_hyp vloc_nonrel_dft_hyp.f90)
target_link_libraries(vloc_nonrel_dft_hyp dftatom)
add_test(vloc_nonrel_dft_hyp ${PROJECT_BINARY_DIR}/vloc_nonrel_dft_hyp)
