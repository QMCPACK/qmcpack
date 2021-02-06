
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)

# E_scf_QP=-150.95914112726956
#  Reference from VMC runs
  LIST(APPEND C4_SD_AE_QP_SHORT "totenergy" "-150.960566064824 0.074761412109")
  LIST(APPEND C4_SD_AE_QP_SHORT "variance" "34.741879231108 4.833758192799")
  LIST(APPEND C4_SD_AE_QP_SHORT "kinetic" "149.744668494694 2.753622966294")
  LIST(APPEND C4_SD_AE_QP_SHORT "potential" "-300.705234559488 2.761778303546")
  LIST(APPEND C4_SD_AE_QP_SHORT "eeenergy" "106.589716452578 0.093871118179")
  LIST(APPEND C4_SD_AE_QP_SHORT "ionion" "56.960598986000 0.000000000000")

  QMC_RUN_AND_CHECK(short-C4_sd-ae-qp-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae"
                    qmc_short-SD-C4_AE_Mol_QP
                    qmc_short-SD-C4_AE_Mol_QP.in.xml 
                    1 16
                    TRUE
                    0 C4_SD_AE_QP_SHORT # VMC
                    )



  LIST(APPEND C4_SD_AE_QP_LONG "totenergy" "-150.960566064824 0.007790118606")
  LIST(APPEND C4_SD_AE_QP_LONG "variance" "34.741879231108 0.990741221699")
  LIST(APPEND C4_SD_AE_QP_LONG "kinetic" "149.744668494694 0.331776551096")
  LIST(APPEND C4_SD_AE_QP_LONG "potential" "-300.705234559488 0.330983245767")
  LIST(APPEND C4_SD_AE_QP_LONG "eeenergy" "106.589716452578 0.011179248515")
  LIST(APPEND C4_SD_AE_QP_LONG "ionion" "56.960598986000 0.000000000000")

  QMC_RUN_AND_CHECK(long-C4_sd-ae-qp-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae"
                    qmc_long-SD-C4_AE_Mol_QP
                    qmc_long-SD-C4_AE_Mol_QP.in.xml 
                    1 16
                    TRUE
                    0 C4_SD_AE_QP_LONG # VMC
                    )

#Nb_Det = 211
#Ground:
#E               =   -151.05936895837109     
#Vnn             =    56.960598986007938     
#Ven             =   -464.01242694605406     
#Vee             =    106.35022537160214     
#T               =    149.64423363007290     
 

#Excited:
#E               =   -150.89066382584013     
#Vnn             =    56.960598986007938     
#Ven             =   -463.18948426388619     
#Vee             =    106.03958987195335     
#T               =    149.30063158008477     


  LIST(APPEND C4_MD_AE_EXCITED_QP_SHORT "totenergy" "-150.892309364193 0.062212972577")
  LIST(APPEND C4_MD_AE_EXCITED_QP_SHORT "variance" "34.255149340568 3.877696647688")
  LIST(APPEND C4_MD_AE_EXCITED_QP_SHORT "kinetic" "149.395366601073 2.289442218166")
  LIST(APPEND C4_MD_AE_EXCITED_QP_SHORT "potential" "-300.287675965416 2.302220898617")
  LIST(APPEND C4_MD_AE_EXCITED_QP_SHORT "eeenergy" "106.039594332889 0.113866937600")
  LIST(APPEND C4_MD_AE_EXCITED_QP_SHORT "ionion" "56.960598986000 0.000000000000")


  QMC_RUN_AND_CHECK(short-C4_md-ae-excited-qp-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae"
                    qmc_short-MD-C4_AE_Mol_Excited_QP
                    qmc_short-MD-C4_AE_Mol_Excited_QP.in.xml
                    1 16
                    TRUE
                    0 C4_MD_AE_EXCITED_QP_SHORT # VMC
                    )


  LIST(APPEND C4_MD_AE_GROUND_QP_SHORT "totenergy" "-151.060022705978 0.067531501861")
  LIST(APPEND C4_MD_AE_GROUND_QP_SHORT "variance" "34.274239480694 1.882676106875")
  LIST(APPEND C4_MD_AE_GROUND_QP_SHORT "kinetic" "149.806496130140 2.275152959288")
  LIST(APPEND C4_MD_AE_GROUND_QP_SHORT "potential" "-300.866518836168 2.280577422331")
  LIST(APPEND C4_MD_AE_GROUND_QP_SHORT "eeenergy" "106.350153396650 0.098463197040")
  LIST(APPEND C4_MD_AE_GROUND_QP_SHORT "ionion" "56.960598986000 0.000000000000")


  QMC_RUN_AND_CHECK(short-C4_md-ae-ground-qp-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae"
                    qmc_short-MD-C4_AE_Mol_Ground_QP
                    qmc_short-MD-C4_AE_Mol_Ground_QP.in.xml
                    1 16
                    TRUE
                    0 C4_MD_AE_GROUND_QP_SHORT # VMC
                    )
 
#--- Deterministic
  
 
 
  LIST(APPEND DET_C4_SD_AE_QP_SHORT "kinetic" "54.11526393 0.00001")
  LIST(APPEND DET_C4_SD_AE_QP_SHORT "totenergy" "-152.43823035 0.00001")
  LIST(APPEND DET_C4_SD_AE_QP_SHORT "potential" "-206.55349428 0.00001")
  LIST(APPEND DET_C4_SD_AE_QP_SHORT "eeenergy" "99.57907161 0.00001")
  LIST(APPEND DET_C4_SD_AE_QP_SHORT "ionion" "56.96059899 0.00001")
  LIST(APPEND DET_C4_SD_AE_QP_SHORT "samples" "9 0.0")
 
  QMC_RUN_AND_CHECK(deterministic-C4_sd-ae-qp-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae"
                    det_qmc_short-SD-C4_AE_Mol_QP
                    det_qmc_short-SD-C4_AE_Mol_QP.in.xml
                    1 1
                    TRUE
                    0 DET_C4_SD_AE_QP_SHORT # VMC
                    )
 
 

  LIST(APPEND DET_C4_MD_AE_GROUND_QP_SHORT "kinetic" "70.64698241 0.00001")
  LIST(APPEND DET_C4_MD_AE_GROUND_QP_SHORT "totenergy" "-154.62718031 0.00001")
  LIST(APPEND DET_C4_MD_AE_GROUND_QP_SHORT "potential" "-225.27416271 0.00001")
  LIST(APPEND DET_C4_MD_AE_GROUND_QP_SHORT "eeenergy" "105.93504749 0.00001")
  LIST(APPEND DET_C4_MD_AE_GROUND_QP_SHORT "ionion" "56.96059899 0.00001")
  QMC_RUN_AND_CHECK(deterministic-C4_md-ae-ground-qp-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae"
                    det_qmc_short-MD-C4_AE_Mol_Ground_QP
                    det_qmc_short-MD-C4_AE_Mol_Ground_QP.in.xml
                     1 1
                    TRUE
                    0 DET_C4_MD_AE_GROUND_QP_SHORT # VMC
                   )
 
 
  LIST(APPEND DET_C4_MD_AE_EXICTED_QP_SHORT "kinetic" "-4.57501459 0.00001")
  LIST(APPEND DET_C4_MD_AE_EXCITED_QP_SHORT "totenergy" "-147.01870207 0.00001")
  LIST(APPEND DET_C4_MD_AE_EXCITED_QP_SHORT "potential" "-142.44368748 0.00001")
   LIST(APPEND DET_C4_MD_AE_EXCITED_QP_SHORT "eeenergy" "95.13749252 0.00001")
  LIST(APPEND DET_C4_MD_AE_EXCITED_QP_SHORT "ionion" "56.96059899 0.00001")
  QMC_RUN_AND_CHECK(deterministic-C4_md-ae-excited-qp-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae"
                    det_qmc_short-MD-C4_AE_Mol_Excited_QP
                    det_qmc_short-MD-C4_AE_Mol_Excited_QP.in.xml
                    1 1
                    TRUE
                    0 DET_C4_MD_AE_EXCITED_QP_SHORT # VMC
                    )





 ELSE()
    MESSAGE_VERBOSE("Skipping C4_ae tests because Pure real WF are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE_VERBOSE("Skipping C4_ae tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

