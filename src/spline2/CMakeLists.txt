# add apps XYZ.cpp, e.g., qmc_particles.cpp
SET(ESTEST bspline_smp bspline_tiledThreading bspline_tiledThreading_compact)

FOREACH(p ${ESTEST})
  #use MultiBsplineStd.hpp : SoA + SIMD/QPX
  ADD_EXECUTABLE( ${p}  tests/${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} ${QMC_UTIL_LIBS}  qmcutil ${MPI_LIBRARY})

  #use MultiBsplineStd2.hpp : Unrolled + SoA + SIMD
  ADD_EXECUTABLE(${p}_std2  tests/${p}.cpp)
  set_target_properties(${p}_std2 PROPERTIES COMPILE_FLAGS "-DUSE_EINSPLINE_UNROLLED")
  TARGET_LINK_LIBRARIES(${p}_std2 ${QMC_UTIL_LIBS}  qmcutil ${MPI_LIBRARY})

  #use MultiBsplineStd4.hpp : AoS + SIMD outer loop
  ADD_EXECUTABLE(${p}_std4  tests/${p}.cpp)
  set_target_properties(${p}_std4 PROPERTIES COMPILE_FLAGS "-DUSE_EINSPLINE_STD4")
  TARGET_LINK_LIBRARIES(${p}_std4 ${QMC_UTIL_LIBS}  qmcutil ${MPI_LIBRARY})

  #use MultiBsplineBase.hpp : AoS
  ADD_EXECUTABLE(${p}_base  tests/${p}.cpp)
  set_target_properties(${p}_base PROPERTIES COMPILE_FLAGS "-DUSE_EINSPLINE_BASE")
  TARGET_LINK_LIBRARIES(${p}_base ${QMC_UTIL_LIBS}  qmcutil ${MPI_LIBRARY})

  #use MultiBsplineBase.hpp : AoS + SIMD
  ADD_EXECUTABLE(${p}_base_simd  tests/${p}.cpp)
  set_target_properties(${p}_base_simd PROPERTIES COMPILE_FLAGS "-DUSE_EINSPLINE_BASE -DENABLE_EINSPLINE_SIMD")
  TARGET_LINK_LIBRARIES(${p}_base_simd ${QMC_UTIL_LIBS}  qmcutil ${MPI_LIBRARY})

ENDFOREACH(p ${ESTEST})

