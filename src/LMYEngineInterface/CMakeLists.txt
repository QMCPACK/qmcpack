# find the LMYEngine library
#FIND_LIBRARY(LMYENGINE_LIB lmyengine)
MESSAGE("-- LMYENGINE_LIB=${LMYENGINE_LIB}")

set(LMY_LIBS "${LMYENGINE_LIB}")

#ADD_DEFINITIONS(-DHAVE_LMY_ENGINE)
# find the LMYEngine header file add it to the includes
#FIND_PATH(LMYENGINE_INCLUDE_DIR formic/fqmc/engine/var_dependencies.h)
#MESSAGE("-- LMYENGINE_INCLUDE_DIR=${LMYENGINE_INCLUDE_DIR}")
#INCLUDE_DIRECTORIES(${LMYENGINE_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_library(lmy_engine_interface LMYE_QMCCostFunctionOMP.cpp
                                 LMYE_QMCFixedSampleLinearOptimize.cpp
           )

set_property(TARGET lmy_engine_interface PROPERTY CXX_STANDARD 11)
set_property(TARGET lmy_engine_interface PROPERTY CXX_STANDARD_REQUIRED ON)

# find the atlas libraries that lmy depends on
#find_library(LMY_ATLAS_LIB_ATLAS atlas)
#find_library(LMY_ATLAS_LIB_CBLAS cblas)
#find_library(LMY_ATLAS_LIB_F77BLAS f77blas)
#find_library(LMY_ATLAS_LIB_LAPACK lapack)
#if(LMY_ATLAS_LIB_ATLAS AND LMY_ATLAS_LIB_CBLAS AND LMY_ATLAS_LIB_F77BLAS AND LMY_ATLAS_LIB_LAPACK)
#  message(STATUS "Found ATLAS for LMY")
#  set(LMY_LIBS "${LMY_LIBS};${LMY_ATLAS_LIB_LAPACK};${LMY_ATLAS_LIB_F77BLAS};${LMY_ATLAS_LIB_CBLAS};${LMY_ATLAS_LIB_ATLAS}")
#else(LMY_ATLAS_LIB_ATLAS AND LMY_ATLAS_LIB_CBLAS AND LMY_ATLAS_LIB_F77BLAS AND LMY_ATLAS_LIB_LAPACK)
#endif(LMY_ATLAS_LIB_ATLAS AND LMY_ATLAS_LIB_CBLAS AND LMY_ATLAS_LIB_F77BLAS AND LMY_ATLAS_LIB_LAPACK)

#if(CMAKE_COMPILER_IS_GNUCXX)
#  set(LMY_LIBS ${LMY_LIBS};m;gfortran)
#endif(CMAKE_COMPILER_IS_GNUCXX)

target_link_libraries(lmy_engine_interface ${LMY_LIBS})

#install(
#        TARGETS
#          lmy_engine_interface
#        RUNTIME DESTINATION
#          bin
#        LIBRARY DESTINATION
#          lib
#        ARCHIVE DESTINATION
#          lib
#       )
