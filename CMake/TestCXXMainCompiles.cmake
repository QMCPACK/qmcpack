# Check that the configured compiler works on a C++ main function
# Note: whitespaces not allowed in STAGE_NAME
function(TestCXXMainCompiles STAGE_NAME)
  if(STAGE_NAME MATCHES " ")
    message(FATAL_ERROR "TestCXXMainCompiles whitespaces not allowed in the stage name. The given value is '${STAGE_NAME}'.")
  endif()
  set(TEST_CXX_COMPILE_MAIN_DIR ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp)
  file(WRITE ${TEST_CXX_COMPILE_MAIN_DIR}/try_cxx_main.cpp "int main(){}")
  set(TEST_RESULT_VAR_NAME TEST_RESULT_${STAGE_NAME})
  try_compile(
    ${TEST_RESULT_VAR_NAME}
    ${TEST_CXX_COMPILE_MAIN_DIR}
    SOURCES ${TEST_CXX_COMPILE_MAIN_DIR}/try_cxx_main.cpp
    OUTPUT_VARIABLE COMPILE_OUTPUT)
  if(NOT ${TEST_RESULT_VAR_NAME})
    message(FATAL_ERROR "Failed in compiling a main() function in stage ${STAGE_NAME}. Output:\n${COMPILE_OUTPUT}")
  endif()
endfunction()
