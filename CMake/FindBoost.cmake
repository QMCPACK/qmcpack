#
# this module look for boost (http://www.boost.org) support
# it will define the following values
#
# BOOST_INCLUDE_DIR = where boost/boost.h can be found
#
# May want to define this but seldom required
# BOOST_LIBRARY = where boost library can be found (reserved)
#
IF(Boost_INCLUDE_DIRS)
  FIND_PATH(BOOST_INCLUDE_DIR boost/config.hpp  ${Boost_INCLUDE_DIRS})
ELSE(Boost_INCLUDE_DIRS)
  SET(TRIAL_PATHS
    $ENV{BOOST_HOME}
    /usr/apps/include
    /usr/include
    /opt/include
    /usr/local/include
  )
  FIND_PATH(BOOST_INCLUDE_DIR boost/config.hpp ${TRIAL_PATHS})
ENDIF(Boost_INCLUDE_DIRS)

IF(BOOST_INCLUDE_DIR)
  SET(FOUND_BOOST 1 CACHE BOOL "Found boost library")
ELSE(BOOST_INCLUDE_DIR)
  SET(FOUND_BOOST 0 CACHE BOOL "Found boost library")
ENDIF(BOOST_INCLUDE_DIR)

MARK_AS_ADVANCED(
   BOOST_INCLUDE_DIR
   FOUND_BOOST
)
