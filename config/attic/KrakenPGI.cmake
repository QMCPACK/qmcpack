# tool chain for Jaguar, Cray XT5  using PrgEnv-gnu
# Fri Feb 20 15:50:10 EST 2009 by Jeongnim Kim
# module list
# Currently Loaded Modulefiles:
# 1) modules/3.1.6                           13) Base-opts/2.1.41HD
# 2) DefApps                                 14) subversion/1.5.0
# 3) torque/2.3.2-snap.200807092141          15) gcc/4.2.0.quadcore
# 4) moab/5.2.4                              16) fftw/3.1.1
# 5) xtpe-quadcore                           17) xt-libsci/10.3.1
# 6) MySQL/5.0.45                            18) xt-mpt/3.1.0
# 7) xt-service/2.1.41HD                     19) xt-pe/2.1.41HD
# 8) xt-libc/2.1.41HD                        20) xt-asyncpe/2.0
# 9) xt-os/2.1.41HD                          21) PrgEnv-gnu/2.1.41HD
# 10) xt-boot/2.1.41HD                        22) cmake/2.6.1
# 11) xt-lustre-ss/2.1.UP00_ORNL.nic52_1.6.5  23) acml/4.1.0
# 12) xtpe-target-cnl                         24) hdf5/1.6.8

SET(CMAKE_SYSTEM_PROCESSOR "XT5")
SET(QMC_ENV "CrayXTEnv" CACHE STRING "Setting envirnoments for Cray XT5")

SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)

set(CMAKE_C_COMPILER  /opt/cray/xt-asyncpe/2.0/bin/cc)
set(CMAKE_CXX_COMPILER  /opt/cray/xt-asyncpe/2.0/bin/CC)

  set(CMAKE_FIND_ROOT_PATH
      /opt/fftw/3.2.0
      /nics/b/home/jnkim2/svnwork/p3dfft/fftw3.2.0
     )

SET(CMAKE_SHARED_LIBRARY_C_FLAGS "")            # -pic 
SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "")       # -shared
SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")         # +s, flag for exe link to use shared lib
SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG "")       # -rpath
SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG_SEP "")   # : or empty

SET(CMAKE_LINK_LIBRARY_SUFFIX "")
SET(CMAKE_STATIC_LIBRARY_PREFIX "lib")
SET(CMAKE_STATIC_LIBRARY_SUFFIX ".a")
SET(CMAKE_SHARED_LIBRARY_PREFIX "lib")          # lib
SET(CMAKE_SHARED_LIBRARY_SUFFIX ".a")           # .a
SET(CMAKE_EXECUTABLE_SUFFIX "")          # .exe
SET(CMAKE_DL_LIBS "" )

SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

#SET(CMAKE_C_FLAGS "-restrict -fastsse -Minline=levels:10")
#SET(CMAKE_CXX_FLAGS "-restrict -fastsse -Minline=levels:10 --no_exceptions -DBOOST_NO_EXCEPTIONS")
SET(CMAKE_C_FLAGS " -c99 -fastsse -Minline=levels:10 -mp=nonuma")
SET(CMAKE_CXX_FLAGS "--restrict  -fastsse -Minline=levels:10 --no_exceptions -DBOOST_NO_EXCEPTIONS -mp=nonuma")
SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS " ")    # -shared
SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS " ")  # +s, flag for exe link to use shared lib

INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)
