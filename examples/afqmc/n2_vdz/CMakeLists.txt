

INCLUDE( "${qmcpack_SOURCE_DIR}/CMake/macros.cmake" )
LIST(APPEND AFQMC_N2_SCALARS "Ebound" "-109.26 0.02")

QMC_RUN_AND_CHECK(short-afqmc-N2_vdz
                  "${CMAKE_SOURCE_DIR}/examples/afqmc/n2_vdz"
                  N2
                  n2.xml
                  4 1
                  AFQMC_N2_SCALARS
                  0 # VMC
                  TRUE)
