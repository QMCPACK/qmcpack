#
# Solids
#

add_subdirectory("LiH_solid_1x1x1_pp")
add_subdirectory("bccH_3x3x3_ae")
add_subdirectory("bccH_2x2x2_ae")
add_subdirectory("bccH_1x1x1_ae")
add_subdirectory("diamondC_1x1x1_pp")
add_subdirectory("diamondC_2x1x1_pp")
add_subdirectory("diamondC_1x1x1-Gaussian_pp")
add_subdirectory("diamondC_1x1x1-Gaussian_pp_MSD")
add_subdirectory("diamondC_2x1x1-Gaussian_pp_MSD")
add_subdirectory("diamondC_2x1x1-Gaussian_pp")
add_subdirectory("diamondC_2x1x1-Gaussian_pp_kpts")
add_subdirectory("hcpBe_1x1x1_pp")
add_subdirectory("monoO_1x1x1_pp")
add_subdirectory("NiO_a4_e48_pp")
add_subdirectory("grapheneC_1x1_pp")
add_subdirectory("InSe_a4")

if(QMC_COMPLEX)
  add_subdirectory("bccMo_2x1x1_SOREP")
  add_subdirectory("diamondC_1x1x1-Gaussian_pp_Tw_cplx")
  add_subdirectory("diamondC_2x1x1-Gaussian_pp_cplx")
  add_subdirectory("diamondC_2x1x1-Gaussian_pp_cplx_MSD")
  add_subdirectory("diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited")
  add_subdirectory("diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD")
  add_subdirectory("Al-1x1x1-Gaussian_pp_ShiftedTwist")
endif()
if(BUILD_AFQMC
   AND QMC_COMPLEX
   AND NOT ENABLE_CUDA)
  add_subdirectory("diamondC_afqmc_1x1x1_complex")
endif()
