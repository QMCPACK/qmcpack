PROJECT(miniapps)

if(QMC_BUILD_LEVEL GREATER 4)
# add apps XYZ.cpp, e.g., qmc_particles.cpp
SET(ESTEST einspline_smp einspline_spo qmc_particles moveonsphere twobody ptclset)
SET(ESTEST diff_distancetables diff_jeeI diff_j2 diff_j1 einspline_spo einspline_spo_nested determinant restart)

FOREACH(p ${ESTEST})
  ADD_EXECUTABLE( ${p}  ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})
ENDFOREACH(p ${ESTEST})

ADD_EXECUTABLE(diff_wfs diff_wfs.cpp AoSWaveFunction.cpp SoAWaveFunction.cpp)
TARGET_LINK_LIBRARIES(diff_wfs qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

ADD_EXECUTABLE(miniqmc_gr miniqmc.cpp AoSWaveFunction.cpp SoAWaveFunction.cpp)
set_target_properties (miniqmc_gr PROPERTIES COMPILE_FLAGS "-DNX=48 -DNY=48 -DNZ=60")
TARGET_LINK_LIBRARIES(miniqmc_gr qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

ADD_EXECUTABLE(miniqmc_nio miniqmc.cpp AoSWaveFunction.cpp SoAWaveFunction.cpp)
set_target_properties (miniqmc_nio PROPERTIES COMPILE_FLAGS "-DUSE_NIO -DNX=37 -DNY=37 -DNZ=37")
TARGET_LINK_LIBRARIES(miniqmc_nio qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

#Use a denser grid for graphite
#ADD_EXECUTABLE(miniqmc_gr2 miniqmc.cpp AoSWaveFunction.cpp SoAWaveFunction.cpp)
#set_target_properties (miniqmc_gr2 PROPERTIES COMPILE_FLAGS "-DNX=64 -DNY=64 -DNZ=80")
#TARGET_LINK_LIBRARIES(miniqmc_gr2 qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

endif()

#SET(boost_test exchange_walker)
#FOREACH(p ${boost_test})
#  ADD_EXECUTABLE( ${p}  ${p}.cpp)
#  TARGET_LINK_LIBRARIES(${p} qmcbase qmcutil ${QMC_UTIL_LIBS} boost ${MPI_LIBRARY})
#ENDFOREACH(p ${boost_test})

