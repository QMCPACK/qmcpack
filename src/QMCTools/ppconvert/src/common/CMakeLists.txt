add_library(common EXCLUDE_FROM_ALL "")

target_sources(common
  PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/GKIntegration.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/CoulombPot.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/NLPP.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/DFTAtom.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/ScreenedPot.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/Fitting.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/RadialWF.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/SplinePot.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/CubicSplineCommon.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/GeneralPot.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/Potential.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/IO.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/MatrixOps.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/Grid.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/PotentialBase.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/IOASCII.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/Functionals.cc"
  )

# get_property( COMMON_COMPILE_DEFS TARGET common PROPERTY COMPILE_DEFINITIONS )
# set(COMMON_COMPILE_DEFS ${COMMON_COMPILE_DEFS} $<$<NOT:$<CXX_COMPILER_ID:Clang>>:__GNUC__=4>)
# MESSAGE( "$<CXX_COMPILER_ID>" )

  target_compile_definitions(common
    PRIVATE
    $<$<CXX_COMPILER_ID:Intel>:__ICC=1>
  )

# "$<$<CXX_COMPILER_ID:Intel>:__ICC>"
#   "$<$<CXX_COMPILER_ID:XL>:__IBM>"
#   "$<$<CXX_COMPILER_ID:PGI>:__PGI>"
#   "$<$<CXX_COMPILER_ID:GNU>:__GNUC__>"

get_property( COMMON_COMPILE_DEFS TARGET common PROPERTY COMPILE_DEFINITIONS )
MESSAGE( "common compile defs ${COMMON_COMPILE_DEFS}" )
target_include_directories(common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/blitz
  ${CMAKE_CURRENT_SOURCE_DIR}/blitz/array
)

