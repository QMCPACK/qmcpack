# Heavily inspired by the CTestCoverageCollectGCOV CMake module.
# Functions for handling coverage data
#  - create gcov files
#  - compare multiple runs
#  - create tar file for CDash upload

# Generate gcov files from gcda and gcno files
# Create the data.json file cdash expects
function(GENERATE_GCOV BINARY_DIR OUTPUT_DIR GCOV_OPTIONS SOURCE_DIR)
  file(MAKE_DIRECTORY ${OUTPUT_DIR})

  file(GLOB_RECURSE GCDA_FILES "${BINARY_DIR}/*.gcda")

  set(GCOV_CMD_OPTIONS "-b;-p")
  if(GCOV_OPTIONS STREQUAL "USE_LONG_FILE_NAMES")
    set(GCOV_CMD_OPTIONS "${GCOV_CMD_OPTIONS};-l;-s;${SOURCE_DIR}")
  endif()
  message("GCOV_CMD_OPTIONS = ${GCOV_CMD_OPTIONS}")

  foreach(GCDA_FILE ${GCDA_FILES})
    execute_process(
      COMMAND gcov ${GCOV_CMD_OPTIONS} ${GCDA_FILE}
      WORKING_DIRECTORY ${OUTPUT_DIR}
      OUTPUT_VARIABLE out)
  endforeach()

  file(
    WRITE ${OUTPUT_DIR}/data.json
    "{
        \"Source\":\"${CTEST_SOURCE_DIRECTORY}\",
        \"Binary\":\"${CTEST_BINARY_DIRECTORY}\"
}")

endfunction()

# Remove unwanted gcov files (files in /usr, unit tests, files with coverage only in static initializers, etc.)
function(FILTER_GCOV GCOV_DIR)
  execute_process(COMMAND ${CTEST_SOURCE_DIRECTORY}/tests/coverage/compare_gcov.py --action process --base-dir
                          ${GCOV_DIR})

endfunction()

# Use after running gcov with the -l (--long-file-names) option to merge all the
#  gcov files from the input directory into one gcov file for each source file in
#  the output directory.
function(MERGE_GCOV INPUT_DIR OUTPUT_DIR SOURCE_DIR)
  file(MAKE_DIRECTORY ${OUTPUT_DIR})

  execute_process(COMMAND python ${CTEST_SOURCE_DIRECTORY}/tests/coverage/compare_gcov.py --action merge --base-dir
                          ${INPUT_DIR} --output-dir ${OUTPUT_DIR} --prefix ${SOURCE_DIR})

  file(COPY ${INPUT_DIR}/data.json DESTINATION ${OUTPUT_DIR})

endfunction()

# Create tar file of gcov files
function(CREATE_GCOV_TAR BINARY_DIRECTORY OUTPUT_DIR)
  execute_process(COMMAND tar cfj gcov_${OUTPUT_DIR}.tar "--mtime=1970-01-01 0:0:0 UTC" ${OUTPUT_DIR}
                  WORKING_DIRECTORY ${BINARY_DIRECTORY})
endfunction()

# Clear the coverage data files in preparation for another run
function(CLEAR_GCDA BINARY_DIRECTORY)
  file(GLOB_RECURSE GCDA_FILES "${BINARY_DIRECTORY}/*.gcda")
  foreach(GCDA_FILE ${GCDA_FILES})
    file(REMOVE ${GCDA_FILE})
  endforeach()
endfunction()

# Compare two coverage runs
function(COMPARE_GCOV BASE_DIR UNIT_DIR OUTPUT_DIR REL_OUTPUT_DIR)
  file(MAKE_DIRECTORY ${OUTPUT_DIR})

  execute_process(COMMAND python ${CTEST_SOURCE_DIRECTORY}/tests/coverage/compare_gcov.py --action compare --base-dir
                          ${BASE_DIR} --unit-dir ${UNIT_DIR} --output-dir ${OUTPUT_DIR})

  file(
    WRITE ${OUTPUT_DIR}/data.json
    "{
        \"Source\":\"${CTEST_SOURCE_DIRECTORY}\",
        \"Binary\":\"${CTEST_BINARY_DIRECTORY}\"
}")

  #FILE(RELATIVE_PATH REL_OUTPUT_DIR ${CTEST_BINARY_DIRECTORY} ${OUTPUT_DIR})
  ##MESSAGE("*** Relative output dir = ${REL_OUTPUT_DIR}")

  execute_process(COMMAND tar cfj gcov.tar "--mtime=1970-01-01 0:0:0 UTC" ${REL_OUTPUT_DIR}
                  WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY})

endfunction()
