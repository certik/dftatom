include_directories(${PROJECT_BINARY_DIR}/src)

project(atom_U)

add_executable(F_atom_U F_atom_U.f90)
target_link_libraries(F_atom_U dftatom)
add_test(F_atom_U ${PROJECT_BINARY_DIR}/F_atom_U)

if(WITH_C_INTERFACE)
    set(CMAKE_C_FLAGS "-std=c99 -Wall -Wextra -fPIC")
    add_executable(C_atom_U C_atom_U.c)
    target_link_libraries(C_atom_U dftatom)
    add_test(C_atom_U ${PROJECT_BINARY_DIR}/C_atom_U)
endif()
