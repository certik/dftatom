#!/usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

if [[ "${TEST_MANUAL}" == "yes" ]]; then
    make -f Makefile.manual test
else
    ctest -E conv --output-on-failure
    if [[ "${WITH_PYTHON}" == "yes" ]]; then
        pytest
        rm -r dftatom
        mkdir xx
        cd xx
        pytest --pyargs dftatom
        python ../examples/atom_B.py
        python ../examples/atom_U.py
    fi
fi
