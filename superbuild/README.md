This uses CMake's [externalproject_add](https://cmake.org/cmake/help/latest/module/ExternalProject.html) to build, test
and install qmcpack with both `-DQMC_COMPLEX=0` and `-DQMC_COMPLEX=1`. Any arguments passed on the command line are passed 
to both builds. To build, change `cmake [options] <path-to-source>` to `cmake [options] <path-to-source>/superbuild`.

Ctest is configured with two tests, one for the real build and another for the complex build. Passing arguments proved 
to be challenging. By default unit tests are run. To change arguments set the environment variable `TEST_ARGUMENTS`. We 
recommend adding the `--verbose` argument to view individual tests.

TESTS_ARGUMENTS="-R short -LE unstable" ctest --verbose

The install target copies all executables from the real build along with `qmcpack_complex` and `qmcpack_complex.settings`
from the complex build.

