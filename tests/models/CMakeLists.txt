# The errors on the total energy are numerical noise, not statistical errors
list(APPEND SHO_SCALARS "totenergy" "0.45  0.00000004")
list(APPEND SHO_SCALARS "kinetic" "0.225 0.0028")
list(APPEND SHO_SCALARS "potential" "0.225 0.0028")
list(APPEND SHO_SCALARS "latdev" "1.0   0.013")

# The errors on the total energy are statistical errors
list(APPEND GEP_SCALARS "totenergy" "0.45  1e-6")
list(APPEND GEP_SCALARS "kinetic" "0.225 0.0030")
list(APPEND GEP_SCALARS "potential" "0.225 0.0030")
list(APPEND GEP_SCALARS "latdev" "1.0   0.013")

list(APPEND DET_SHO_VMC_SCALARS "totenergy" "0.4500000000 0.000001")
list(APPEND DET_SHO_VMC_SCALARS "kinetic" "0.3142864882 0.000001")
list(APPEND DET_SHO_VMC_SCALARS "potential" "0.1357135118 0.000001")
list(APPEND DET_SHO_VMC_SCALARS "latdev" "0.6031711637 0.000001")

list(APPEND DET_SHO_ALLP_VMC_SCALARS "totenergy" "0.4500000000 0.000001")
list(APPEND DET_SHO_ALLP_VMC_SCALARS "kinetic" "0.3142864882 0.000001")
list(APPEND DET_SHO_ALLP_VMC_SCALARS "potential" "0.1357135118 0.000001")
list(APPEND DET_SHO_ALLP_VMC_SCALARS "latdev" "0.6031711637 0.000001")

list(APPEND DET_SHO_NODRIFT_VMC_SCALARS "totenergy" "0.4500000000 0.000001")
list(APPEND DET_SHO_NODRIFT_VMC_SCALARS "kinetic" "0.3545249251 0.000001")
list(APPEND DET_SHO_NODRIFT_VMC_SCALARS "potential" "0.0954750749 0.000001")
list(APPEND DET_SHO_NODRIFT_VMC_SCALARS "latdev" "0.4243336664 0.000001")

list(APPEND DET_SHO_NODRIFT_VMCBATCH_SCALARS "totenergy" "0.4500000000 0.000001")
list(APPEND DET_SHO_NODRIFT_VMCBATCH_SCALARS "kinetic" "0.3545249251 0.000001")
list(APPEND DET_SHO_NODRIFT_VMCBATCH_SCALARS "potential" "0.0954750749 0.000001")
list(APPEND DET_SHO_NODRIFT_VMCBATCH_SCALARS "latdev" "0.4243336664 0.000001")

list(APPEND DET_SHO_VMCBATCH_SCALARS "totenergy" "0.4500000000 0.000001")
list(APPEND DET_SHO_VMCBATCH_SCALARS "kinetic" "0.3142864882 0.000001")
list(APPEND DET_SHO_VMCBATCH_SCALARS "potential" "0.1357135118 0.000001")
list(APPEND DET_SHO_VMCBATCH_SCALARS "latdev" "0.6031711637 0.000001")

list(APPEND DET_GEP_VMC_SCALARS "totenergy" "0.4500000396 0.000001")
list(APPEND DET_GEP_VMC_SCALARS "kinetic" "0.3142864882 0.000001")
list(APPEND DET_GEP_VMC_SCALARS "potential" "0.1357135514 0.000001")
list(APPEND DET_GEP_VMC_SCALARS "latdev" "0.6031711637 0.000001")

list(APPEND DET_GEP_VMCBATCH_SCALARS "totenergy" "0.4500000666 0.000001")
list(APPEND DET_GEP_VMCBATCH_SCALARS "kinetic" "0.3142864882 0.000001")
list(APPEND DET_GEP_VMCBATCH_SCALARS "potential" "0.1357135118 0.000001")
list(APPEND DET_GEP_VMCBATCH_SCALARS "latdev" "0.6031711637 0.000001")

qmc_run_and_check(
  short-sho-vmc
  "${qmcpack_SOURCE_DIR}/tests/models/sho"
  vmc_sho
  vmc_sho.xml
  1
  16
  TRUE
  0
  SHO_SCALARS # VMC
)

qmc_run_and_check(
  short-sho-vmc-allp
  "${qmcpack_SOURCE_DIR}/tests/models/sho"
  vmc_sho_allp
  vmc_sho_allp.xml
  1
  16
  TRUE
  0
  SHO_SCALARS # VMC
)

qmc_run_and_check(
  short-sho-vmc-nodrift
  "${qmcpack_SOURCE_DIR}/tests/models/sho"
  vmc_sho_nodrift
  vmc_sho_nodrift.xml
  1
  16
  TRUE
  0
  SHO_SCALARS # VMC
)


qmc_run_and_check(
  deterministic-sho-vmc
  "${qmcpack_SOURCE_DIR}/tests/models/sho"
  det_vmc_sho
  det_vmc_sho.xml
  1
  1
  TRUE
  0
  DET_SHO_VMC_SCALARS
)

qmc_run_and_check(
  deterministic-sho-vmc-allp
  "${qmcpack_SOURCE_DIR}/tests/models/sho"
  det_vmc_sho_allp
  det_vmc_sho_allp.xml
  1
  1
  TRUE
  0
  DET_SHO_ALLP_VMC_SCALARS
)


qmc_run_and_check(
  deterministic-sho-vmc-nodrift
  "${qmcpack_SOURCE_DIR}/tests/models/sho"
  det_vmc_sho_nodrift
  det_vmc_sho_nodrift.xml
  1
  1
  TRUE
  0
  DET_SHO_NODRIFT_VMC_SCALARS
)

qmc_run_and_check(
  short-sho-vmcbatch
  "${qmcpack_SOURCE_DIR}/tests/models/sho"
  vmcbatch_sho
  vmcbatch_sho.xml
  1
  16
  TRUE
  0
  SHO_SCALARS # VMC
)

qmc_run_and_check(
  short-sho-vmcbatch-nodrift 
  "${qmcpack_SOURCE_DIR}/tests/models/sho"
  vmcbatch_sho_nodrift
  vmcbatch_sho_nodrift.xml
  1
  16
  TRUE
  0
  SHO_SCALARS # VMC
)

qmc_run_and_check(
  deterministic-sho-vmcbatch
  "${qmcpack_SOURCE_DIR}/tests/models/sho"
  det_vmcbatch_sho
  det_vmcbatch_sho.xml
  1
  1
  TRUE
  0
  DET_SHO_VMCBATCH_SCALARS
)

qmc_run_and_check(
  deterministic-sho-vmcbatch-nodrift
  "${qmcpack_SOURCE_DIR}/tests/models/sho"
  det_vmcbatch_sho_nodrift
  det_vmcbatch_sho_nodrift.xml
  1
  1
  TRUE
  0
  DET_SHO_NODRIFT_VMCBATCH_SCALARS
)

simple_run_and_check(deterministic-sho-grad_lap "${qmcpack_SOURCE_DIR}/tests/models/sho" grad_lap.xml 1 1 check_grad_lap.py)

qmc_run_and_check(
  short-grid_ext_pot-vmc
  "${qmcpack_SOURCE_DIR}/tests/models/grid_external_potential"
  vmc_grid_ext_pot
  vmc_grid_ext_pot.xml
  1
  16
  TRUE
  0
  GEP_SCALARS # VMC
)

qmc_run_and_check(
  deterministic-grid_ext_pot-vmc
  "${qmcpack_SOURCE_DIR}/tests/models/grid_external_potential"
  det_vmc_grid_ext_pot
  det_vmc_grid_ext_pot.xml
  1
  1
  TRUE
  0
  DET_GEP_VMC_SCALARS
)

qmc_run_and_check(
  short-grid_ext_pot-vmcbatch
  "${qmcpack_SOURCE_DIR}/tests/models/grid_external_potential"
  vmcbatch_grid_ext_pot
  vmcbatch_grid_ext_pot.xml
  1
  16
  TRUE
  0
  GEP_SCALARS # VMC
)

qmc_run_and_check(
  deterministic-grid_ext_pot-vmcbatch
  "${qmcpack_SOURCE_DIR}/tests/models/grid_external_potential"
  det_vmcbatch_grid_ext_pot
  det_vmcbatch_grid_ext_pot.xml
  1
  1
  TRUE
  0
  DET_GEP_VMCBATCH_SCALARS
)
