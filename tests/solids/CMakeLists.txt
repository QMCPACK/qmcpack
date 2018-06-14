

#
# Solids
#

SUBDIRS("LiH_solid_1x1x1_pp")
SUBDIRS("bccH_3x3x3_ae")
SUBDIRS("bccH_2x2x2_ae")
SUBDIRS("bccH_1x1x1_ae")
SUBDIRS("diamondC_1x1x1_pp")
SUBDIRS("diamondC_2x1x1_pp")
SUBDIRS("diamondC_1x1x1-Gaussian_pp")
SUBDIRS("diamondC_2x1x1-Gaussian_pp")
SUBDIRS("hcpBe_1x1x1_pp")
SUBDIRS("monoO_1x1x1_pp")
SUBDIRS("NiO_a4_e48_pp")
IF (BUILD_AFQMC AND NOT QMC_COMPLEX)
   SUBDIRS("diamondC_afqmc")
   SUBDIRS("NiO_afqmc")
ENDIF()

