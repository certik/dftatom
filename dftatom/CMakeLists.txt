add_subdirectory(lib)

set(PY_PATH ${PYTHON_INSTALL_PATH}/dftatom)
install(FILES __init__.py
    DESTINATION ${PY_PATH}
    )
install(FILES tests/test_dftatom.py
    DESTINATION ${PY_PATH}/tests
    )
