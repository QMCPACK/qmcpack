# Test example

IF(NOT QMC_CUDA AND NOT QMC_COMPLEX)

  LIST(APPEND HE_SIMPLE_VMC_SCALARS "totenergy" "-2.83 .006") # total energy

  QMC_RUN_AND_CHECK(example_He_simple
                    "${CMAKE_SOURCE_DIR}/examples/molecules/He"
                    He 
                    he_simple.xml
                    1 1
                    HE_SIMPLE_VMC_SCALARS
                    0
                    TRUE)

  LIST(APPEND HE_SIMPLE_DMC_SCALARS "totenergy" "-2.902 .002") # total energy

  QMC_RUN_AND_CHECK(example_He_simple_dmc
                    "${CMAKE_SOURCE_DIR}/examples/molecules/He"
                    He 
                    he_simple_dmc.xml
                    1 1
                    HE_SIMPLE_DMC_SCALARS
                    1 # series for DMC data
                    TRUE)

  LIST(APPEND HE_SIMPLE_OPT_SCALARS "totenergy" "-2.88 .004") # total energy

  QMC_RUN_AND_CHECK(example_He_simple_opt
                    "${CMAKE_SOURCE_DIR}/examples/molecules/He"
                    He 
                    he_simple_opt.xml
                    1 1
                    HE_SIMPLE_OPT_SCALARS
                    10 # final optimized series
                    TRUE)

  LIST(APPEND HE_FROM_GAMESS_SCALARS "totenergy" "-2.872 .003") # total energy

  QMC_RUN_AND_CHECK(example_He_from_gamess
                    "${CMAKE_SOURCE_DIR}/examples/molecules/He"
                    He 
                    he_from_gamess.xml
                    1 1
                    HE_FROM_GAMESS_SCALARS
                    0
                    TRUE)

# Not included for now because it fails spectacularly a small percentage of
# the time (2/100 runs).  Values at each optimization stage seem to be okay,
# but the final VMC is far off.
#  LIST(APPEND HE_BSPLINE_OPT_SCALARS "totenergy" "-2.879 .003") # total energy
#
#  QMC_RUN_AND_CHECK(example_He_bspline_opt
#                    "${CMAKE_SOURCE_DIR}/examples/molecules/He"
#                    He 
#                    he_bspline_jastrow.xml
#                    1 1
#                    HE_BSPLINE_OPT_SCALARS
#                    10 # series for VMC on optimized value
#                    TRUE)
ENDIF()
