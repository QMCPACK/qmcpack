PROJECT(pp)

#CONFIGURE_FILE(${hf_SOURCE_DIR}/../Configuration.h.in
#              ${hf_SOURCE_DIR}/../Configuration.h)

SET(PPSRCS 
  ../Utilities/OhmmsInform.cpp 
  ../OhmmsApp/ProjectData.cpp 
  ../Platforms/sysutil.cpp
  ../AtomicHF/Clebsch_Gordan.cpp
  ../AtomicHF/RadialPotential.cpp
  ../PseudoGen/PseudoGen.cpp
)

LINK_LIBRARIES(${FORTRAN_LIBS})
ADD_EXECUTABLE(pp ${PPSRCS})
TARGET_LINK_LIBRARIES(pp qmcbase)
