include_directories(${PROJECT_BINARY_DIR}/src)

project(enrichment)
add_executable(enrichment enrichment.f90)
target_link_libraries(enrichment dftatom)
