IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)

#Nb_Det = 1
# -150.95531020041986

#  Reference from VMC runs
  LIST(APPEND C4_SD_AE_NO_SHORT "totenergy" "-150.95310200 0.07367645")

  QMC_RUN_AND_CHECK(short-C4_sd-ae-no-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
                    qmc_short-SD-C4_AE_Mol_NO
                    qmc_short-SD-C4_AE_Mol_NO.in.xml 
                    1 16
                    TRUE
                    0 C4_SD_AE_NO_SHORT # VMC
                    )

#Nb_Det = 160
#Ground:
#E               =   -151.08986108933223
#Vnn             =    56.960598881874070
#Ven             =   -464.83401896934589
#Vee             =    106.77210949007733
#T               =    150.01144950806227

#Excited:
#E               =   -150.93124043303607
#Vnn             =    56.960598881874070
#Ven             =   -465.62023598802489
#Vee             =    107.09690302643341
#T               =    150.63149364668132

  LIST(APPEND C4_MD_AE_GROUND_NO_SHORT "kinetic" "150.01144950806227 2.94012241")
  LIST(APPEND C4_MD_AE_GROUND_NO_SHORT "totenergy" "-151.08986108933223 0.06871005")
  LIST(APPEND C4_MD_AE_GROUND_NO_SHORT "potential" "-301.10131059739449 2.92475134")
  LIST(APPEND C4_MD_AE_GROUND_NO_SHORT "eeenergy" "106.77210949007733 0.08571242")
  LIST(APPEND C4_MD_AE_GROUND_NO_SHORT "ionion" "56.960598881874070 0.0003")
  LIST(APPEND C4_MD_AE_GROUND_NO_SHORT "samples" "128000 0.0")

  QMC_RUN_AND_CHECK(short-C4_md-ae-ground-no-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
                    qmc_short-MD-C4_AE_Mol_Ground_NO
                    qmc_short-MD-C4_AE_Mol_Ground_NO.in.xml
                    1 16
                    TRUE
                    0 C4_MD_AE_GROUND_NO_SHORT # VMC
                    )

  LIST(APPEND C4_MD_AE_EXCITED_NO_SHORT "kinetic" "150.63149364668132 3.22648038")
  LIST(APPEND C4_MD_AE_EXCITED_NO_SHORT "totenergy" "-150.93124043303607 0.08723064")
  LIST(APPEND C4_MD_AE_EXCITED_NO_SHORT "potential" "-301.56273407971741 3.24423157")
  LIST(APPEND C4_MD_AE_EXCITED_NO_SHORT "eeenergy" "107.09690302643341 0.09393176")
  LIST(APPEND C4_MD_AE_EXCITED_NO_SHORT "ionion" "56.960598881874070 0.0003")
  LIST(APPEND C4_MD_AE_EXCITED_NO_SHORT "samples" "128000 0.0")

  QMC_RUN_AND_CHECK(short-C4_md-ae-excited-no-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
                    qmc_short-MD-C4_AE_Mol_Excited_NO
                    qmc_short-MD-C4_AE_Mol_Excited_NO.in.xml
                    1 16
                    TRUE
                    0 C4_MD_AE_EXCITED_NO_SHORT # VMC
                    )

#Nb_Det = 1
# -150.95914112726956

  LIST(APPEND C4_SD_AE_QP_SHORT "totenergy" "-150.95914112726956 0.05956643")

  QMC_RUN_AND_CHECK(short-C4_sd-ae-qp-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
                    qmc_short-SD-C4_AE_Mol_QP
                    qmc_short-SD-C4_AE_Mol_QP.in.xml
                    1 16
                    TRUE
                    0 C4_SD_AE_QP_SHORT # VMC
                    )


#Nb_Det = 211
#Ground:
#E               =   -151.05736895837109
#Vnn             =    56.960598986007938
#Ven             =   -464.01242694605406
#Vee             =    106.35022537160214
#T               =    149.64423363007290


#Excited:
#E               =   -150.88866382584013
#Vnn             =    56.960598986007938
#Ven             =   -463.18948426388619
#Vee             =    106.03958987195335
#T               =    149.30063158008477

#  Reference from VMC runs
  LIST(APPEND C4_MD_AE_GROUND_QP_SHORT "kinetic" "149.64423363007290 1.94055367")
  LIST(APPEND C4_MD_AE_GROUND_QP_SHORT "totenergy" "-151.05736895837109 0.06818342")
  LIST(APPEND C4_MD_AE_GROUND_QP_SHORT "potential" "-300.701602588443982 1.93622152")
  LIST(APPEND C4_MD_AE_GROUND_QP_SHORT "eeenergy" "106.35022537160214 0.10457402")
  LIST(APPEND C4_MD_AE_GROUND_QP_SHORT "ionion" "56.960598986007938 0.0003")

  QMC_RUN_AND_CHECK(short-C4_md-ae-ground-qp-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
                    qmc_short-MD-C4_AE_Mol_Ground_QP
                    qmc_short-MD-C4_AE_Mol_Ground_QP.in.xml
                    1 16
                    TRUE
                    0 C4_MD_AE_GROUND_QP_SHORT # VMC
                   )


  LIST(APPEND C4_MD_AE_EXICTED_QP_SHORT "kinetic" "149.30063158008477 4.63130094")
  LIST(APPEND C4_MD_AE_EXCITED_QP_SHORT "totenergy" "-150.88866382584013 0.06554183")
  LIST(APPEND C4_MD_AE_EXCITED_QP_SHORT "potential" "-300.189295405924902 4.60708030")
  LIST(APPEND C4_MD_AE_EXCITED_QP_SHORT "eeenergy" "106.03958987195335 0.09946428")
  LIST(APPEND C4_MD_AE_EXCITED_QP_SHORT "ionion" "56.960598986007938 0.0003")

  QMC_RUN_AND_CHECK(short-C4_md-ae-excited-qp-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
                    qmc_short-MD-C4_AE_Mol_Excited_QP
                    qmc_short-MD-C4_AE_Mol_Excited_QP.in.xml
                    1 16
                    TRUE
                    0 C4_MD_AE_EXCITED_QP_SHORT # VMC
                   )


#Deterministic tests

  LIST(APPEND DET_C4_SD_AE_NO_SHORT "kinetic" "40.88568798 0.00001") 
  LIST(APPEND DET_C4_SD_AE_NO_SHORT "totenergy" "-146.21104723 0.00001") 
  LIST(APPEND DET_C4_SD_AE_NO_SHORT "potential" "-187.09673521 0.00001")
  LIST(APPEND DET_C4_SD_AE_NO_SHORT "eeenergy" "99.68606745 0.00001")
  LIST(APPEND DET_C4_SD_AE_NO_SHORT "ionion" "56.96059899 0.00001")
  LIST(APPEND DET_C4_SD_AE_NO_SHORT "samples" "9 0.0") 
  QMC_RUN_AND_CHECK(deterministic-C4_sd-ae-no-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
                    det_qmc_short-SD-C4_AE_Mol_NO
                    det_qmc_short-SD-C4_AE_Mol_NO.in.xml
                    1 1
                    TRUE
                    0 DET_C4_SD_AE_NO_SHORT # VMC
                    )


  LIST(APPEND DET_C4_MD_AE_GROUND_NO_SHORT "kinetic" "81.37875523 0.00001") 
  LIST(APPEND DET_C4_MD_AE_GROUND_NO_SHORT "totenergy" "-147.35211603 0.00001") 
  LIST(APPEND DET_C4_MD_AE_GROUND_NO_SHORT "potential" "-228.73087126 0.00001")
  LIST(APPEND DET_C4_MD_AE_GROUND_NO_SHORT "eeenergy" "94.66385499 0.00001") 
  LIST(APPEND DET_C4_MD_AE_GROUND_NO_SHORT "ionion" "56.96059899 0.00001")
  LIST(APPEND DET_C4_MD_AE_GROUND_NO_SHORT "samples" "9 0.0") 
  QMC_RUN_AND_CHECK(deterministic-C4_md-ae-ground-no-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
                    det_qmc_short-MD-C4_AE_Mol_Ground_NO
                    det_qmc_short-MD-C4_AE_Mol_Ground_NO.in.xml
                    1 1
                    TRUE
                    0 DET_C4_MD_AE_GROUND_NO_SHORT # VMC
                    )


  LIST(APPEND DET_C4_MD_AE_EXCITED_NO_SHORT "kinetic" "120.73710114 0.00001")
  LIST(APPEND DET_C4_MD_AE_EXCITED_NO_SHORT "totenergy" "-139.80104570 0.00001")
  LIST(APPEND DET_C4_MD_AE_EXCITED_NO_SHORT "potential" "-260.53814684 0.00001")
  LIST(APPEND DET_C4_MD_AE_EXCITED_NO_SHORT "eeenergy" "114.94712080 0.00001")
  LIST(APPEND DET_C4_MD_AE_EXCITED_NO_SHORT "ionion" "56.96059899 0.00001")
  LIST(APPEND DET_C4_MD_AE_EXCITED_NO_SHORT "samples" "9 0.0")
  QMC_RUN_AND_CHECK(deterministic-C4_md-ae-excited-no-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
                    det_qmc_short-MD-C4_AE_Mol_Excited_NO
                    det_qmc_short-MD-C4_AE_Mol_Excited_NO.in.xml
                    1 1
                    TRUE
                    0 DET_C4_MD_AE_EXCITED_NO_SHORT # VMC
                    )



  LIST(APPEND DET_C4_SD_AE_QP_SHORT "kinetic" "54.11526393 0.00001")
  LIST(APPEND DET_C4_SD_AE_QP_SHORT "totenergy" "-152.43823035 0.00001")
  LIST(APPEND DET_C4_SD_AE_QP_SHORT "potential" "-206.55349428 0.00001")
  LIST(APPEND DET_C4_SD_AE_QP_SHORT "eeenergy" "99.57907161 0.00001")
  LIST(APPEND DET_C4_SD_AE_QP_SHORT "ionion" "56.96059899 0.00001")
  LIST(APPEND DET_C4_SD_AE_QP_SHORT "samples" "9 0.0")

  QMC_RUN_AND_CHECK(deterministic-C4_sd-ae-qp-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
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
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
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
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C4_ae_MSD_excited"
                    det_qmc_short-MD-C4_AE_Mol_Excited_QP
                    det_qmc_short-MD-C4_AE_Mol_Excited_QP.in.xml
                    1 1
                    TRUE
                    0 DET_C4_MD_AE_EXCITED_QP_SHORT # VMC
                   )



 ELSE()
    MESSAGE_VERBOSE("Skipping C2_MSD_Excited tests because Pure real WF are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE_VERBOSE("Skipping C2_MSD_Excited tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()
