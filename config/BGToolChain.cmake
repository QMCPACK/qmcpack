# the name of the target operating system
#SET(CMAKE_SYSTEM_NAME BlueGeneP)
SET(BGP 1 CACHE BOOL "On BlueGeneP")
SET(Linux 0)

# set the compiler
#set(CMAKE_C_COMPILER  /opt/ibmcmp/vacpp/bg/9.0/bin/bgxlc_r)
#set(CMAKE_CXX_COMPILER  /opt/ibmcmp/vacpp/bg/9.0/bin/bgxlC_r)
set(CMAKE_C_COMPILER  /bgsys/drivers/ppcfloor/comm/bin/mpixlc_r)
set(CMAKE_CXX_COMPILER  /bgsys/drivers/ppcfloor/comm/bin/mpixlcxx_r)

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
    /home/straka/build
    /home/straka/jk2/ohmmscore/CMake
)

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_SHARED_LINKER_FLAGS " ")

set(AIX_ARCH "450")
#SET(AIX_ARCH_FLAGS "-qarch=${AIX_ARCH}  -qtune=${AIX_ARCH}  -qcache=${AIX_ARCH} -qsmp=omp")
SET(AIX_ARCH_FLAGS "-qarch=${AIX_ARCH}  -qtune=${AIX_ARCH} -qsmp=omp")
#SET(AIX_CXX_COMMON_FLAGS " -qkeyword=restrict -qstrict -qhot -qrtti=dyna -qtemplaterecompile -qnoeh -qsuppress=1540-1090 ")
SET(AIX_CXX_COMMON_FLAGS " -qkeyword=restrict -qstrict -qhot -qtemplaterecompile -qnoeh -qsuppress=1540-1090:1540-1088 ")
SET(AIX_OPT_FLAGS "-O3 -Q -qmaxmem=-1 -qipa=inline -qinline -qlargepage -qprefetch")
#SET(AIX_CXX_OPT_FLAGS "-O3 -Q -qlargepage -qprefetch")
#SET(AIX_CXX_FLAGS "-O3 -Q -qlargepage -qprefetch")

SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_CXX_COMMON_FLAGS} ${AIX_OPT_FLAGS}")
SET(CMAKE_C_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS}")

#SET(CMAKE_CXX_CREATE_STATIC_LIBRARY
#       "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS> "
#       "<CMAKE_RANLIB> <TARGET> ")

SET(CMAKE_C_LINK_EXECUTABLE
"<CMAKE_C_COMPILER> <FLAGS> <CMAKE_C_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> "
"  -o <TARGET> <LINK_LIBRARIES> -Wl,-lc ")

SET(CMAKE_CXX_LINK_EXECUTABLE
"<CMAKE_CXX_COMPILER> <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> "
" -o <TARGET> <LINK_LIBRARIES> -Wl,-lstdc++")

