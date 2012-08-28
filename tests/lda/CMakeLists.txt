include_directories(${PROJECT_BINARY_DIR}/src)

project(nist)

add_executable(conv_lda conv_lda.f90)
target_link_libraries(conv_lda dftatom)

add_test(conv_lda ${PROJECT_BINARY_DIR}/conv_lda)
