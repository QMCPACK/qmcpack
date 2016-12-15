# Ensure that examples run

SUBDIRS(molecules/H2O)
SUBDIRS(molecules/He)
IF (WITH_AFQMC)
  SUBDIRS(afqmc/n2_vdz)
ENDIF()
