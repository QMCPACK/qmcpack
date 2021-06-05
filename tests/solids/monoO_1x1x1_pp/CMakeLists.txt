list(APPEND MONO_O_SCALARS "totenergy" "-31.776596 0.006")
list(APPEND MONO_O_SCALARS "samples" "64000 0.0")
#  LIST(APPEND MONO_O_SCALARS "flux" "0.0 0.5")

qmc_run_and_check(
  short-monoO_1x1x1_pp-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
  qmc_short
  qmc_short.in.xml
  1
  16
  TRUE
  0
  MONO_O_SCALARS # VMC
)

list(APPEND MONO_O_J3_SCALARS "totenergy" "-31.913838 0.0025")
list(APPEND MONO_O_J3_SCALARS "samples" "204800 0.0")
qmc_run_and_check(
  short-monoO_1x1x1_pp-vmc_sdj3
  "${qmcpack_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
  qmc_j3_short
  qmc_j3_short.in.xml
  1
  16
  TRUE
  0
  MONO_O_J3_SCALARS # VMC
)

#
# Long test
#

list(APPEND LONG_MONO_O_SCALARS "totenergy" "-31.776596 0.0002")
list(APPEND LONG_MONO_O_SCALARS "samples" "61440000 0.0")
list(APPEND LONG_MONO_O_SCALARS "flux" "0.0 0.05")

qmc_run_and_check(
  long-monoO_1x1x1_pp-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
  qmc_long
  qmc_long.in.xml
  1
  16
  TRUE
  0
  LONG_MONO_O_SCALARS # VMC
)

list(APPEND LONG_MONO_O_J3_SCALARS "totenergy" "-31.913838 0.0008")
list(APPEND LONG_MONO_O_J3_SCALARS "samples" "2048000 0.0")
list(APPEND LONG_MONO_O_J3_SCALARS "flux" "0.0 0.27")

qmc_run_and_check(
  long-monoO_1x1x1_pp-vmc_sdj3
  "${qmcpack_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
  qmc_j3_long
  qmc_j3_long.in.xml
  1
  16
  TRUE
  0
  LONG_MONO_O_J3_SCALARS # VMC
)
