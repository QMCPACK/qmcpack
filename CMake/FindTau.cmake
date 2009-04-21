#
# this module look for tau (http://www.boost.org) support
# it will define the following values
#
# TAU_INCLUDE_DIR = where TAU.h can be found
#
# May want to define this but seldom required
# TAU_LIBRARY = where boost library can be found (reserved)
#
SET(TRIAL_PATHS
    ${CMAKE_FIND_ROOT_PATH}
    /usr/apps/include
    /usr/include
    /opt/include
    /usr/local/include
   )

IF($ENV{TAUROOT} MATCHES "tau")
  SET(TAU_HOME $ENV{TAUROOT})
ENDIF($ENV{TAUROOT} MATCHES "tau")

IF($ENV{TAU_HOME} MATCHES "tau")
  SET(TAU_HOME $ENV{TAU_HOME})
ENDIF($ENV{TAU_HOME} MATCHES "tau")

SET(TRIAL_PATHS ${TRIAL_PATHS} ${TAU_HOME}/include)

FIND_PATH(TAU_INCLUDE_DIR TAU.h ${TRIAL_PATHS})

IF(TAU_INCLUDE_DIR)
  SET(FOUND_TAU 1 CACHE BOOL "Found tau library")
ELSE(TAU_INCLUDE_DIR)
  SET(FOUND_TAU 0 CACHE BOOL "Found tau library")
ENDIF(TAU_INCLUDE_DIR)

MARK_AS_ADVANCED(
  TAU_INCLUDE_DIR
  FOUND_TAU
)
