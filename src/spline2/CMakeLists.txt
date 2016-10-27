# add apps XYZ.cpp, e.g., qmc_particles.cpp
SET(ESTEST bspline_smp bspline_tiledThreading bspline_tiledThreading_compact)

FOREACH(p ${ESTEST})
  #use MultiBsplineStd.hpp : SoA + SIMD/QPX
  ADD_EXECUTABLE( ${p}  tests/${p}.cpp)
  target_link_libraries(${p} qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

  #use MultiBsplineStd2.hpp : Unrolled + SoA + SIMD
  ADD_EXECUTABLE(${p}_std2  tests/${p}.cpp)
  set_target_properties(${p}_std2 PROPERTIES COMPILE_FLAGS "-DUSE_EINSPLINE_UNROLLED")
  target_link_libraries(${p}_std2 qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})
  
  #use MultiBsplineStd4.hpp : AoS + SIMD outer loop
  ADD_EXECUTABLE(${p}_std4  tests/${p}.cpp)
  set_target_properties(${p}_std4 PROPERTIES COMPILE_FLAGS "-DUSE_EINSPLINE_STD4")
  target_link_libraries(${p}_std4 qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})
  
  #use MultiBsplineBase.hpp : AoS
  ADD_EXECUTABLE(${p}_base  tests/${p}.cpp)
  set_target_properties(${p}_base PROPERTIES COMPILE_FLAGS "-DUSE_EINSPLINE_BASE")
  target_link_libraries(${p}_base qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})
  
  #use MultiBsplineBase.hpp : AoS + SIMD
  ADD_EXECUTABLE(${p}_base_simd  tests/${p}.cpp)
  set_target_properties(${p}_base_simd PROPERTIES COMPILE_FLAGS "-DUSE_EINSPLINE_BASE -DENABLE_EINSPLINE_SIMD")
  target_link_libraries(${p}_base_simd qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

ENDFOREACH(p ${ESTEST})

