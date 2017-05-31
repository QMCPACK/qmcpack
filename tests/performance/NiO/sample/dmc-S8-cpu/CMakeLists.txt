
SET(PERF_TEST NiO-cpu-S8-dmc)
SET(TEST_DIR dmc-S8-cpu)
SET(INPUT_FILE NiO-fcc-S8-dmc.xml)
SET(H5_FILE NiO-fcc-supertwist111-supershift000-S8.h5)

ADD_NIO_TEST(${PERF_TEST} ${TEST_DIR} ${INPUT_FILE} ${H5_FILE})


# performance tests
#SET(PERF_TEST NiO-cpu-S8-dmc)
#
#SET (WDIR "${qmcpack_BINARY_DIR}/tests/performance/NiO/dmc-S8-cpu")
#SET (SDIR "${CMAKE_SOURCE_DIR}/tests/performance/NiO/dmc-S8-cpu")
#
#SET(QMC_APP "${qmcpack_BINARY_DIR}/bin/qmcpack")
#SET(PERF_ARGS "--enable-timers=fine")
#
#SET(INPUT_FILE NiO-fcc-S8-dmc.xml)
#
#SET(H5_FILE NiO-fcc-supertwist111-supershift000-S8.h5)
#
#SET (PPDIR "${CMAKE_SOURCE_DIR}/tests/performance/NiO")
#
#EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory ${SDIR} ${WDIR})
#
#FILE(GLOB PP "${PPDIR}/*.xml")
#MESSAGE("pp files = ${PP}")
#FOREACH(F IN LISTS PP)
#  FILE(RELATIVE_PATH R ${PPDIR} ${F})
#  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink "${F}" ${WDIR}/${R})
#ENDFOREACH()
#
#EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink "${QMC_DATA_NIO_H5}/${H5_FILE}" "${WDIR}/${H5_FILE}")
#
# Need the INPUT_FILE if grabbing the whole thing externally.  Should this
# have a particular form (.xml? .in.xml?)
#ADD_TEST(NAME ${PERF_TEST} COMMAND ${QMC_APP} ${PERF_ARGS} ${INPUT_FILE})
#SET_TESTS_PROPERTIES(${PERF_TEST} PROPERTIES LABELS "performance")
#SET_TESTS_PROPERTIES(${PERF_TEST} PROPERTIES WORKING_DIRECTORY "${WDIR}")
