set(SRC
    dftatom_wrapper.cpp
    )

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
include_directories(${dftatom_BINARY_DIR}/src)

cython_add_module_pyx(dftatom_wrapper)
add_python_library(dftatom_wrapper ${SRC})
target_link_libraries(dftatom_wrapper dftatom)

add_custom_target(cython
    COMMAND cython dftatom_wrapper.pyx
    )

set(PY_PATH ${PYTHON_INSTALL_PATH}/dftatom/lib)
install(TARGETS dftatom_wrapper
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
