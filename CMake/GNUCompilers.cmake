#GNU compilers
IF(CMAKE_COMPILER_IS_GNUCXX) 
  ADD_DEFINITIONS(-Drestrict=__restrict__ -DADD_)
  SET(CMAKE_CXX_FLAGS "-O6 -ftemplate-depth-60 -Drestrict=__restrict__ -fstrict-aliasing -funroll-all-loops   -finline-limit=1000 -ffast-math -Wno-deprecated")
#  SET(CMAKE_CXX_FLAGS "-g -ftemplate-depth-60 -Drestrict=__restrict__ -fstrict-aliasing -Wno-deprecated")
  SET(FORTRAN_LIBS "-lg2c")
  SET(F77 g77)
  SET(F77FLAGS  -funroll-loops -O3)
ENDIF(CMAKE_COMPILER_IS_GNUCXX) 
  
IF(APPLE)
  INCLUDE_DIRECTORIES(/sw/include)
ENDIF(APPLE)
