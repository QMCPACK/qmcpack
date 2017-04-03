
# Heavily inspired by the CTestCoverageCollectGCOV CMake module.
# Functions for handling coverage data
#  - create gcov files
#  - compare multiple runs
#  - create tar file for CDash upload

# Generate gcov files from gcda and gcno files
# Create the data.json file cdash expects
FUNCTION(GENERATE_GCOV BINARY_DIR OUTPUT_DIR OUTPUT_BASE)
  FILE(MAKE_DIRECTORY ${OUTPUT_DIR})

  FILE(GLOB_RECURSE GCDA_FILES "${BINARY_DIR}/*.gcda")

  FOREACH(GCDA_FILE ${GCDA_FILES})
    EXECUTE_PROCESS(COMMAND gcov -b -p ${GCDA_FILE} WORKING_DIRECTORY ${OUTPUT_DIR} OUTPUT_VARIABLE out)
  ENDFOREACH()

  FILE(WRITE ${OUTPUT_DIR}/data.json
      "{
        \"Source\":\"${CTEST_SOURCE_DIRECTORY}\",
        \"Binary\":\"${CTEST_BINARY_DIRECTORY}\"
}")

ENDFUNCTION()


# Create tar file of gcov files
FUNCTION(CREATE_GCOV_TAR BINARY_DIRECTORY OUTPUT_DIR)
  EXECUTE_PROCESS(COMMAND tar cfj gcov_${OUTPUT_DIR}.tar
                  "--mtime=1970-01-01 0:0:0 UTC"
                  ${OUTPUT_DIR}
                  WORKING_DIRECTORY ${BINARY_DIRECTORY})
ENDFUNCTION()


# Clear the coverage data files in preparation for another run
FUNCTION(CLEAR_GCDA BINARY_DIRECTORY)
  FILE(GLOB_RECURSE GCDA_FILES "${BINARY_DIRECTORY}/*.gcda")
  FOREACH(GCDA_FILE ${GCDA_FILES})
    FILE(REMOVE ${GCDA_FILE})
  ENDFOREACH()
ENDFUNCTION()


# Compare two coverage runs
FUNCTION(COMPARE_GCOV BASE_DIR UNIT_DIR OUTPUT_DIR REL_OUTPUT_DIR)
  FILE(MAKE_DIRECTORY ${OUTPUT_DIR})

  EXECUTE_PROCESS(COMMAND python ${CTEST_SOURCE_DIRECTORY}/tests/coverage/compare_gcov.py --action compare --base-dir ${BASE_DIR} --unit-dir ${UNIT_DIR} --output-dir ${OUTPUT_DIR})


  FILE(WRITE ${OUTPUT_DIR}/data.json
      "{
        \"Source\":\"${CTEST_SOURCE_DIRECTORY}\",
        \"Binary\":\"${CTEST_BINARY_DIRECTORY}\"
}")

  #FILE(RELATIVE_PATH REL_OUTPUT_DIR ${CTEST_BINARY_DIRECTORY} ${OUTPUT_DIR})
  ##MESSAGE("*** Relative output dir = ${REL_OUTPUT_DIR}")

  EXECUTE_PROCESS(COMMAND tar cfj gcov.tar
                  "--mtime=1970-01-01 0:0:0 UTC"
                  ${REL_OUTPUT_DIR}
                  WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY})

ENDFUNCTION()

