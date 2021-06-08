list(APPEND HCP_BE_SCALARS "totenergy" "-1.481656 0.0005")
list(APPEND HCP_BE_SCALARS "samples" "416000 0.0")
#  LIST(APPEND HCP_BE_SCALARS "flux" "0.0 0.003")

qmc_run_and_check(
  short-hcpBe_1x1x1_pp-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/hcpBe_1x1x1_pp"
  qmc_short
  qmc_short.in.xml
  1
  16
  TRUE
  0
  HCP_BE_SCALARS # VMC
)

#
# Long test
#

list(APPEND LONG_HCP_BE_SCALARS "totenergy" "-1.481656 0.000017")
list(APPEND LONG_HCP_BE_SCALARS "samples" "399360000 0.0")
list(APPEND LONG_HCP_BE_SCALARS "flux" "0.0 0.0001")

qmc_run_and_check(
  long-hcpBe_1x1x1_pp-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/hcpBe_1x1x1_pp"
  qmc_long
  qmc_long.in.xml
  1
  16
  TRUE
  0
  LONG_HCP_BE_SCALARS # VMC
)
