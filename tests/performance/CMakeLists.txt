#############################################################
# Add performance tests to ctest
#############################################################

IF (NOT QMC_DATA)
  MESSAGE("QMC_DATA is not set. Performance tests will be skipped.")
ELSE()
  MESSAGE("Adding performance tests for QMCPACK")
ENDIF()

# includes 
add_subdirectory(NiO)
add_subdirectory(C-graphite)
add_subdirectory(C-molecule)

