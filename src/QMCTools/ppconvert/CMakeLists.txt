#####################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0)
IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)
###################################################

subdirs(src)
