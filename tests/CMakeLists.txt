SUBDIRS(converter)
SUBDIRS(estimator)
SUBDIRS(models)
SUBDIRS(io)
SUBDIRS(performance)

IF(QE_BIN)
  MESSAGE(STATUS "QE binary path ${QE_BIN}")
  SUBDIRS(pw2qmcpack)
ENDIF()

message("Adding system tests for QMCPACK")

IF (NOT TEST_MAX_PROCS )
SET ( TEST_MAX_PROCS 16 )
message ("TEST_MAX_PROCS was unset. Set to 16")
ENDIF ()

set(QMCAPP_ERR "")
if((NOT ${TEST_MAX_PROCS}) OR (${TEST_MAX_PROCS} STREQUAL ""))
    set(QMCAPP_ERR "${QMCAPP_ERR}  TEST_MAX_PROCS not set: '${TEST_MAX_PROCS}'\n")
endif()

IF(ENABLE_SOA)
  SET(SOA_SUCCESS TRUE)
ELSE()
  SET(SOA_SUCCESS FALSE)
ENDIF()

IF(MIXED_PRECISION)
  SET(MP_SUCCESS FALSE)
ELSE()
  SET(MP_SUCCESS TRUE)
ENDIF()

#
#
if(NOT ${QMCAPP_ERR} STREQUAL "")
    message("${QMCAPP_ERR}  skipping tests")
else()
    INCLUDE( "${qmcpack_SOURCE_DIR}/CMake/macros.cmake" )
    SUBDIRS("heg/heg_14_gamma")
    SUBDIRS("molecules")
    SUBDIRS("solids")
endif()
