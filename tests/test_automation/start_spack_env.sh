#!/bin/bash

# This shared script fragment start spack and loads in QMCPACKS currently supported versions of compilers
# and libraries

. $SPACK_ROOT/share/spack/setup-env.sh

# This is to deal with the fact that we don't know the relationship between the working directory
# and checkout directory in a relative way.
# this doesn't deal with symlinks
SRC="${BASH_SOURCE[0]}"
SRC_DIR="$( cd -P "$( dirname "$SRC" )" >/dev/null 2>&1 && pwd )"

. ${SRC_DIR}/spack_supported_package_versions.sh
