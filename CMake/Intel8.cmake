#compiler flags for intel 8.x
SET(INTEL_COMPILER 1)
ADD_DEFINITIONS(-DADD_)

#enable Interprocedural (IP) Optimizations
#-ipo_obj force generation of real object files (requires -ipo)
#SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=1 -ansi -ipo -ipo_obj")
#SET(CMAKE_CC_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=1 -ansi -ipo -ipo_obj")
SET(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=1 -ansi")
SET(CMAKE_CC_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=1 -ansi")

#IF(BITS MATCHES 64)
#  SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -ftz")
#ELSE(BITS MATCHES 64)
#  SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -xW")
#ENDIF(BITS MATCHES 64)

IF(OHMMS_OMP)
  SET(ENABLE_OPENMP 1)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
ENDIF(OHMMS_OMP)

#check if ia64 or ia32 and add the appropriate flags 
SET(IA32 i686)
EXEC_PROGRAM(/bin/arch OUTPUT_VARIABLE IA32)
IF(IA32 MATCHES ia64)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -ftz")
ELSE(IA32 MATCHES ia64)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -xW")
ENDIF(IA32 MATCHES ia64)

#ifc -> ifort
SET(FORTRAN_LIBS " -lifcore -lguide")
SET(F77 ifort)
SET(F77OPTFLAGS  -fpp2 -O3)
SET(F77FLAGS ${F77OPTFLAGS})
