

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
SUBDIRS("diamondC_2x1x1-Gaussian_pp_kpts")
SUBDIRS("hcpBe_1x1x1_pp")
SUBDIRS("monoO_1x1x1_pp")
SUBDIRS("NiO_a4_e48_pp")
SUBDIRS("grapheneC_1x1_pp")
IF (QMC_COMPLEX)
   SUBDIRS("diamondC_1x1x1-Gaussian_pp_Tw_cplx")
   SUBDIRS("diamondC_2x1x1-Gaussian_pp_cplx")
   SUBDIRS("diamondC_2x1x1-Gaussian_pp_cplx_MSD")
   SUBDIRS("diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited")
   SUBDIRS("diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD")
   SUBDIRS("Al-1x1x1-Gaussian_pp_ShiftedTwist")
ENDIF()
IF (BUILD_AFQMC AND QMC_COMPLEX AND NOT ENABLE_CUDA)
   SUBDIRS("diamondC_afqmc_1x1x1_complex")
ENDIF()

