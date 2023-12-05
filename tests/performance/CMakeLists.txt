#############################################################
# Add performance tests to ctest
#############################################################

if(NOT QMC_DATA)
  message("QMC_DATA is not set. Performance tests will be skipped.")
else()
  message("Adding performance tests for QMCPACK")
endif()

# includes
add_subdirectory(NiO)
add_subdirectory(C-graphite)
add_subdirectory(C-molecule)
