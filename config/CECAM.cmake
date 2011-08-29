SET(CMAKE_SYSTEM_PROCESSOR "core2")


set(CMAKE_C_COMPILER  gcc)
set(CMAKE_CXX_COMPILER  g++)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(GNU_FLAGS "-fopenmp -O3 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(XT_FLAGS "-msse3 -D_CRAYMPI")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS}")

set(CMAKE_FIND_ROOT_PATH
	/nfs_home/tutoadmin/qmc
	/nfs_home/tutoadmin/qmc/libxml2
   )

SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 0)
SET(HAVE_SSE 1)
SET(HAVE_SSE2 1)
SET(HAVE_SSE3 1)
SET(HAVE_SSSE3 1)
SET(USE_PREFETCH 1)
SET(PREFETCH_AHEAD 12)

link_libraries(-llapack -lblas)
