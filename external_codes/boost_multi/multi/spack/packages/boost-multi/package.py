# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class BoostMulti(Package):
    """Boost.Multi: multidimensional array access to contiguous or regularly
    contiguous memory in C++.

    A header-only modern C++ library that provides access and manipulation of
    data in multidimensional arrays, for synthetic types as well. (Not an
    official Boost library.)"""

    homepage = "https://gitlab.com/correaa/boost-multi"
    git = "https://gitlab.com/correaa/boost-multi.git"
    url = "https://gitlab.com/correaa/boost-multi/-/archive/0.90.0/boost-multi-0.90.0.tar.gz"

    maintainers("correaa")

    license("BSL-1.0")

    version("develop", branch="develop")
    # Run `spack checksum boost-multi` to fill in real hashes for tagged releases:
    version("0.90.0", sha256="0000000000000000000000000000000000000000000000000000000000000000")

    # Header-only: declares the C++ usage for dependents; no compiled build deps.
    depends_on("cxx", type="build")

    # Nothing to compile — just install the headers (like the Conan header-library recipe).
    def install(self, spec, prefix):
        install_tree("include", prefix.include)
