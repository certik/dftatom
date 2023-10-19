cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_SOURCE_DIR}/cmake/UserOverride.cmake)

project(dftatom)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)

# Make sure that CMAKE_BUILD_TYPE is either Debug or Release:
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release
        CACHE STRING "Build type (Debug, Release)" FORCE)
endif ()
if (NOT (CMAKE_BUILD_TYPE STREQUAL "Debug" OR
        CMAKE_BUILD_TYPE STREQUAL "Release"))
    message("${CMAKE_BUILD_TYPE}")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of: Debug, Release (current value: '${CMAKE_BUILD_TYPE}')")
endif ()

set(WITH_PYTHON no
    CACHE BOOL "Build with Python wrappers")
set(WITH_C_INTERFACE no
    CACHE BOOL "Build with C interface")
if (WITH_PYTHON)
    set(WITH_C_INTERFACE yes)
endif ()

set(WITH_LAPACK no
    CACHE BOOL "Build tests that depend on Lapack")
if(WITH_LAPACK)
    find_package(Lapack REQUIRED)
    set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
endif()


if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # gfortran
    # Enable this if you want to check for single/double corruption (and use
    # the Debug build):
    #set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fdefault-real-8")
endif ()


enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
if (WITH_PYTHON)
    # We need to execute the find_package() commands here, so that the
    # variables like ${PYTHON_INSTALL_PATH} are defined here and in all cmake
    # files (as opposed to just subset of cmake files).
    find_package(Python REQUIRED)
    find_package(NumPy REQUIRED)
    find_package(Cython REQUIRED)
    add_subdirectory(dftatom)
endif ()

message("\n")
message("Configuration results")
message("---------------------")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message("Build type: ${CMAKE_BUILD_TYPE}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_DEBUG}")
else ()
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_RELEASE}")
endif ()
message("Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message("LIBS: ${LIBS}")
if (WITH_PYTHON)
    message("Python install path: ${PYTHON_INSTALL_PATH}")
endif ()
message("With C interface: ${WITH_C_INTERFACE}")
message("With Python: ${WITH_PYTHON}")
message("With Lapack: ${WITH_LAPACK}")
if (WITH_LAPACK)
    message("LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
endif(WITH_LAPACK)
