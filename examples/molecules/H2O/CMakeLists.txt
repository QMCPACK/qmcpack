# Test example

IF(NOT QMC_CUDA AND NOT QMC_COMPLEX)

  LIST(APPEND H2O_SCALARS "totenergy" "-17.26138 .006") # total energy

  QMC_RUN_AND_CHECK(example_H2O
                    "${CMAKE_SOURCE_DIR}/examples/molecules/H2O"
                    H2O 
                    simple-H2O.xml
                    1 1
                    TRUE
                    2 H2O_SCALARS # check the DMC data
                    )
ENDIF()
