#####################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0)
IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)
###################################################

###################################################
# Done searching libraries. Start building ppconvert
###################################################

MESSAGE(STATUS "Ready to build ppconvert")

subdirs(src)

# add_library(common EXCLUDE_FROM_ALL "")
# add_executable(ppconvert EXCLUDE_FROM_ALL "")

#include
# include(src/CMakeLists.txt)
