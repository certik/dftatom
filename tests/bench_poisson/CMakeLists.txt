include_directories(${PROJECT_BINARY_DIR}/src)

project(bench_poisson)

add_executable(bench_poisson bench_poisson.f90)
target_link_libraries(bench_poisson dftatom)
add_test(bench_poisson ${PROJECT_BINARY_DIR}/bench_poisson)

add_executable(test_poisson_kernels test_poisson_kernels.f90)
target_link_libraries(test_poisson_kernels dftatom)
add_test(test_poisson_kernels ${PROJECT_BINARY_DIR}/test_poisson_kernels)
