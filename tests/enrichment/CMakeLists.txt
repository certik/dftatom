include_directories(${PROJECT_BINARY_DIR}/src)

project(enrichment)
add_executable(enrichment enrichment.f90)
target_link_libraries(enrichment dftatom)

project(enrichment_ks)
add_executable(enrichment_ks enrichment_ks.f90)
target_link_libraries(enrichment_ks dftatom)
