if(NOT QMC_CUDA)
  if(NOT QMC_MIXED_PRECISION)


  list(APPEND HE_BSPLINE_PARAM jud_0 -0.121914 0.0014) # scalar name, value, error
  list(APPEND HE_BSPLINE_PARAM jud_1 0.069689 0.0015)
  list(APPEND HE_BSPLINE_PARAM jud_2 0.051412 0.00031)
  list(APPEND HE_BSPLINE_PARAM jud_3 0.000812 0.000038)

  qmc_run_and_check_custom_scalar(
    BASE_NAME
    He_param_grad
    BASE_DIR
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_param"
    PREFIX
    He_param_grad.param
    INPUT_FILE
    He_param_grad.xml
    PROCS
    1
    THREADS
    16
    SERIES
    0
    SCALAR_VALUES
    HE_BSPLINE_PARAM)

  # Check the parameter gradient code for the legacy driver
  qmc_run_and_check_custom_scalar(
    BASE_NAME
    He_param_grad_legacy_driver
    BASE_DIR
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_param"
    PREFIX
    He_param_grad_legacy_driver.param
    INPUT_FILE
    He_param_grad_legacy_driver.xml
    PROCS
    1
    THREADS
    16
    SERIES
    0
    SCALAR_VALUES
    HE_BSPLINE_PARAM)


  #  Test loading from variational parameter file
  check_python_reqs(h5py he_param_h5 add_h5_tests)
  if (add_h5_tests)
    if (QMC_COMPLEX)
      set(complex_flag "--complex")
    else()
      set(complex_flag)
    endif()
    set(SDIR "${CMAKE_CURRENT_SOURCE_DIR}")
    # Normally the directory is created by qmc_run_and_check_custom_scalar,
    # but we need to write a file there before that runs
    set(TDIR "${CMAKE_CURRENT_BINARY_DIR}/He_param_grad_load-1-16")
    file(MAKE_DIRECTORY ${TDIR})
    execute_process(COMMAND ${Python3_EXECUTABLE} ${qmcpack_SOURCE_DIR}/tests/molecules/He_param/convert_vp_format.py ${SDIR}/he_vp_opt.txt ${complex_flag} -o ${TDIR}/he_vp_opt.h5)

    list(APPEND HE_BSPLINE_OPT_PARAM jud_0  0.00000124  0.0014) # scalar name, value, error
    list(APPEND HE_BSPLINE_OPT_PARAM jud_1 -0.000273    0.00097)
    list(APPEND HE_BSPLINE_OPT_PARAM jud_2 -0.000181    0.00082)
    list(APPEND HE_BSPLINE_OPT_PARAM jud_3  0.0004463   0.000058)

    qmc_run_and_check_custom_scalar(
      BASE_NAME He_param_grad_load
      BASE_DIR "${qmcpack_SOURCE_DIR}/tests/molecules/He_param"
      PREFIX He_param_grad_load.param
      INPUT_FILE He_param_grad_load.xml
      PROCS 1
      THREADS 16
      SERIES 0
      SCALAR_VALUES HE_BSPLINE_OPT_PARAM)

  endif()


  # Orbital Rotation

  # Two atomic orbitals, no jastrow

  list(APPEND HE_ORB_ROT_PARAM spo-up_orb_rot_0000_0001  -0.218  0.01) # scalar name, value, error
  list(APPEND HE_ORB_ROT_PARAM spo-down_orb_rot_0000_0001  -0.218  0.01) # scalar name, value, error

  qmc_run_and_check_custom_scalar(
    BASE_NAME
    He_orb_rot_param_grad_legacy
    BASE_DIR
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_param"
    PREFIX
    He_orb_rot_param_grad_legacy.param
    INPUT_FILE
    He_orb_rot_param_grad_legacy.xml
    PROCS
    1
    THREADS
    1
    SERIES
    0
    SCALAR_VALUES
    HE_ORB_ROT_PARAM)

  qmc_run_and_check_custom_scalar(
    BASE_NAME
    He_orb_rot_param_grad_legacy
    BASE_DIR
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_param"
    PREFIX
    He_orb_rot_param_grad_legacy.param
    INPUT_FILE
    He_orb_rot_param_grad_legacy.xml
    PROCS
    1
    THREADS
    16
    SERIES
    0
    SCALAR_VALUES
    HE_ORB_ROT_PARAM)


  else()
    message(VERBOSE "Skipping He_param tests because parameter output is not supported by mixed precison build (QMC_MIXED_PRECISION=1)")
  endif()

else()
  message(VERBOSE "Skipping He_param tests because parameter output is not supported by CUDA build (QMC_CUDA=1)")
endif()
