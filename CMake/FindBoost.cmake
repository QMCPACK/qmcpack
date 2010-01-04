#
# this module look for boost (http://www.boost.org) support
# it will define the following values
#
# Boost_INCLUDE_DIR = where boost/boost.h can be found
# Boost_FOUND = boolean
#
# May want to define this but seldom required
# BOOST_LIBRARIES = where boost library can be found (reserved)
#
IF(Boost_INCLUDE_DIRS)
  FIND_PATH(BOOST_INCLUDE_DIR boost/config.hpp  ${Boost_INCLUDE_DIRS})
ELSE(Boost_INCLUDE_DIRS)
  if(EXISTS $ENV{BOOST_HOME}/boost/config.hpp)
    set(Boost_INCLUDE_DIR $ENV{BOOST_HOME})
  else()
    FIND_PATH(Boost_INCLUDE_DIR boost/config.hpp ${QMC_INCLUDE_PATHS})
  endif()
ENDIF(Boost_INCLUDE_DIRS)

SET(Boost_FOUND false)
if(Boost_INCLUDE_DIR)
  SET(Boost_FOUND true)
endif(Boost_INCLUDE_DIR)

MARK_AS_ADVANCED(
   Boost_INCLUDE_DIR
   Boost_FOUND
)
