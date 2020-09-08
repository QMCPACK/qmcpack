#!/bin/bash

# Update packages that track master
# Currently only llvm
# Intended to be run nightly/weekly etc.

echo --- Update Script START `date`

if [ -e `dirname "$0"`/ornl_versions.sh ]; then
    source `dirname "$0"`/ornl_versions.sh
else
    echo Did not find version numbers script ornl_versions.sh
    exit 1
fi

spack uninstall -y llvm@master
spack install llvm@master ^python@${python_version}%gcc@${gcc_vnew}

echo --- Update Script END `date`
