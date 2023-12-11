# The errors on the total energy are numerical noise, not statistical errors
if(QMC_MIXED_PRECISION)
  list(APPEND SHO_SCALARS "totenergy" "0.45  0.00000004")
else()
  list(APPEND SHO_SCALARS "totenergy" "0.45  4.0e-14")
endif()
list(APPEND SHO_SCALARS "kinetic" "0.225 0.0028")
list(APPEND SHO_SCALARS "potential" "0.225 0.0028")
list(APPEND SHO_SCALARS "latdev" "1.0   0.013")

# The errors on the total energy are statistical errors
if(QMC_MIXED_PRECISION)
  list(APPEND GEP_SCALARS "totenergy" "0.45  1e-6")
else()
  list(APPEND GEP_SCALARS "totenergy" "0.45  1e-6")
endif()
list(APPEND GEP_SCALARS "kinetic" "0.225 0.0028")
list(APPEND GEP_SCALARS "potential" "0.225 0.0028")
list(APPEND GEP_SCALARS "latdev" "1.0   0.013")

qmc_run_and_check(
  short-sho-vmc
  "${qmcpack_SOURCE_DIR}/tests/models/sho"
  m5
  vmc_sho.xml
  1
  16
  TRUE
  0
  SHO_SCALARS # VMC
)

simple_run_and_check(short-sho-grad_lap "${qmcpack_SOURCE_DIR}/tests/models/sho" grad_lap.xml 1 1 check_grad_lap.py)

qmc_run_and_check(
  short-grid_ext_pot-vmc
  "${qmcpack_SOURCE_DIR}/tests/models/grid_external_potential"
  grid_ext_pot
  vmc_grid_ext_pot.xml
  1
  16
  TRUE
  0
  GEP_SCALARS # VMC
)
