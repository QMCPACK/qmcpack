PROJECT(hf)

#CONFIGURE_FILE(${hf_SOURCE_DIR}/../Configuration.h.in
#              ${hf_SOURCE_DIR}/../Configuration.h)


SET(HFSRCS 
  ../Utilities/OhmmsInform.cpp 
  ../OhmmsApp/ProjectData.cpp 
  ../Platforms/sysutil.cpp
  ../AtomicHF/Clebsch_Gordan.cpp
  ../AtomicHF/RadialPotential.cpp
  ../AtomicHF/HartreeFock.cpp
)

LINK_LIBRARIES(${FORTRAN_LIBS})
ADD_EXECUTABLE(hf ${HFSRCS})
TARGET_LINK_LIBRARIES(hf qmcbase)
