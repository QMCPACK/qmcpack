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

target_include_directories(common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/blitz
  ${CMAKE_CURRENT_SOURCE_DIR}/blitz/array
)
  
