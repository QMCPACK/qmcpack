
# Heavily inspired by the CTestCoverageCollectGCOV CMake module.
# Functions for handling coverage data
#  - create gcov files
#  - compare multiple runs
#  - create tar file for CDash upload

# Generate gcov files from gcda and gcno files
# Create the data.json file cdash expects
FUNCTION(GENERATE_GCOV BINARY_DIR OUTPUT_DIR GCOV_OPTIONS SOURCE_DIR)
  FILE(MAKE_DIRECTORY ${OUTPUT_DIR})

  FILE(GLOB_RECURSE GCDA_FILES "${BINARY_DIR}/*.gcda")

  SET(GCOV_CMD_OPTIONS "-b;-p")
  IF (GCOV_OPTIONS STREQUAL "USE_LONG_FILE_NAMES")
    SET(GCOV_CMD_OPTIONS "${GCOV_CMD_OPTIONS};-l;-s;${SOURCE_DIR}")
  ENDIF()
  MESSAGE("GCOV_CMD_OPTIONS = ${GCOV_CMD_OPTIONS}")

  FOREACH(GCDA_FILE ${GCDA_FILES})
    EXECUTE_PROCESS(COMMAND gcov ${GCOV_CMD_OPTIONS} ${GCDA_FILE} WORKING_DIRECTORY ${OUTPUT_DIR} OUTPUT_VARIABLE out)
  ENDFOREACH()

  FILE(WRITE ${OUTPUT_DIR}/data.json
      "{
        \"Source\":\"${CTEST_SOURCE_DIRECTORY}\",
        \"Binary\":\"${CTEST_BINARY_DIRECTORY}\"
}")

ENDFUNCTION()

# Remove unwanted gcov files (files in /usr, unit tests, files with coverage only in static initializers, etc.)
FUNCTION(FILTER_GCOV GCOV_DIR)
  EXECUTE_PROCESS(COMMAND python ${CTEST_SOURCE_DIRECTORY}/tests/coverage/compare_gcov.py --action process --base-dir ${GCOV_DIR})

ENDFUNCTION()

# Use after running gcov with the -l (--long-file-names) option to merge all the
#  gcov files from the input directory into one gcov file for each source file in
#  the output directory.
FUNCTION(MERGE_GCOV INPUT_DIR OUTPUT_DIR SOURCE_DIR)
  FILE(MAKE_DIRECTORY ${OUTPUT_DIR})

  EXECUTE_PROCESS(COMMAND python ${CTEST_SOURCE_DIRECTORY}/tests/coverage/compare_gcov.py --action merge --base-dir ${INPUT_DIR} --output-dir ${OUTPUT_DIR} --prefix ${SOURCE_DIR})

  FILE(COPY ${INPUT_DIR}/data.json DESTINATION ${OUTPUT_DIR})

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

