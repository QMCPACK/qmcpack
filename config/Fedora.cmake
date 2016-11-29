
#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpic++ )
set(GNU_OPTS "-fopenmp -Wall -O3 -march=native -DADD_ -DINLINE_ALL=inline -DMAX_CHAR_PER_WALKER=98 -D_TIMER_ -D_LINUX_ -DUSE_MPI ")
set(CMAKE_C_FLAGS "${GNU_OPTS} -std=c99")
set(CMAKE_CXX_FLAGS "${GNU_OPTS} -std=c++11 -Wno-sign-compare -Drestrict=__restrict__")

set(ENABLE_OPENMP 1)
set(HAVE_MPI 1)
set(HAVE_SSE 1)
set(HAVE_SSE2 1)
set(HAVE_SSE3 1)
set(HAVE_SSSE3 1)
set(HAVE_SSE41 1)
set(USE_PREFETCH 1)
set(PREFETCH_AHEAD 10)

link_libraries(rt)
link_libraries(lapack)
link_libraries(blas)

set(TEST_MAX_PROCS 8)
set(QMC_MPI 1)
set(HAVE_MPI 1)
