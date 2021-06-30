if(NOT QMC_CUDA)

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
else()
  message_verbose("Skipping He_param tests because parameter output is not supported by CUDA build (QMC_CUDA=1)")
endif()
