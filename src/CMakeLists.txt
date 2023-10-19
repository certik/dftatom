set(SRC
    types.f90
    constants.f90
    utils.f90
    mesh.f90
    ode1d.f90
    rschroed.f90
    rdirac.f90
    rpoisson.f90
    reigen.f90
    dft.f90
    dft_data.f90
    mixings.f90
    states.f90
    drivers.f90
    energies.f90

    ode1d_other.f90
    rschroed_other.f90
    rdirac_other.f90
    rpoisson_other.f90

    dftatom.f90
)

if(WITH_LAPACK)
    set(SRC ${SRC} lapack.f90 interpolation.f90)
endif()

if(WITH_C_INTERFACE)
    set(SRC ${SRC} c_dftatom.f90)
endif()

add_library(dftatom SHARED ${SRC})
target_include_directories(dftatom PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include/dftatom>
)
target_link_libraries(dftatom ${LIBS})
install(TARGETS dftatom
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib
        )
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include/dftatom)
