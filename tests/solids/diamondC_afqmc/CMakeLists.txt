# add short tests
list(APPEND SHORT_SCALARS "Eloc" "-10.55274 0.00261")
list(APPEND SHORT_HYBRID_SCALARS "ElocEstim" "-10.55149 0.00643")

qmc_run_and_check(
  short-diamondC-afqmc
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc"
  qmc_short
  qmc_short.in.xml
  16
  1
  TRUE
  0
  SHORT_SCALARS # AFQMC
)

# distribute cholesky matrix
qmc_run_and_check(
  short-diamondC-afqmc_nn2
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc"
  qmc_short_nnodes2
  qmc_short_nnodes2.in.xml
  16
  1
  TRUE
  0
  SHORT_SCALARS # AFQMC
)

# hybrid weighting
qmc_run_and_check(
  short-diamondC-afqmc_hyb
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc"
  qmc_short_hybrid
  qmc_short_hybrid.in.xml
  16
  1
  TRUE
  0
  SHORT_HYBRID_SCALARS # AFQMC
)

# hybrid weighting with cholesky distribution
qmc_run_and_check(
  short-diamondC-afqmc_hyb_nn2
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc"
  qmc_short_hybrid_nnodes2
  qmc_short_hybrid_nnodes2.in.xml
  16
  1
  TRUE
  0
  SHORT_HYBRID_SCALARS # AFQMC
)

# include mean field
qmc_run_and_check(
  short-diamondC-afqmc_incmf
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc"
  qmc_short_includeMF
  qmc_short_includeMF.in.xml
  16
  1
  TRUE
  0
  SHORT_SCALARS # AFQMC
)

# include mean field with cholesky distribution
qmc_run_and_check(
  short-diamondC-afqmc_incmf_nn2
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc"
  qmc_short_includeMF_nnodes2
  qmc_short_includeMF_nnodes2.in.xml
  16
  1
  TRUE
  0
  SHORT_SCALARS # AFQMC
)

# add long tests
list(APPEND LONG_SCALARS "Eloc" "-10.55274 0.00086")
list(APPEND LONG_HYBRID_SCALARS "ElocEstim" "-10.55149 0.00212")

qmc_run_and_check(
  long-diamondC-afqmc
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc"
  qmc_long
  qmc_long.in.xml
  16
  1
  TRUE
  0
  LONG_SCALARS # AFQMC
)

qmc_run_and_check(
  long-diamondC-afqmc_hyb
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc"
  qmc_long_hybrid
  qmc_long_hybrid.in.xml
  16
  1
  TRUE
  0
  LONG_HYBRID_SCALARS # AFQMC
)
