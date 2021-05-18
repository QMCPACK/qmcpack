cmake_minimum_required(VERSION 3.14.0)

# Pass CMake arguments to external projects
# https://stackoverflow.com/questions/12021448/how-can-cmake-arguments-be-forwarded-to-externalproject

# MUST be done before call to 'project'
get_cmake_property(vars CACHE_VARIABLES)
foreach(var ${vars})
  get_property(currentHelpString CACHE "${var}" PROPERTY HELPSTRING)
    if("${currentHelpString}" MATCHES "No help, variable specified on the command line." OR "${currentHelpString}" STREQUAL "")
	#message("${var} = [${${var}}]  --  ${currentHelpString}") # uncomment to see the variables being processed
        list(APPEND CL_ARGS "-D${var}=${${var}}")
    endif()
endforeach()

project(qmcpacksuperbuild)

# Don't pass CMAKE_INSTALL_PREFIX to external projects
foreach(item ${CL_ARGS})
  string(FIND ${item} "-DCMAKE_INSTALL_PREFIX=" install_prefix)
  if(${install_prefix} GREATER -1)
    list(REMOVE_ITEM CL_ARGS ${item})
    break()
  endif()
endforeach()

# Don't pass QMC_COMPLEX to external projects
foreach(item ${CL_ARGS})
  string(FIND ${item} "-DQMC_COMPLEX=" complex_option)
  if(${complex_option} GREATER -1)
    message(WARNING "Do not specify QMC_COMPLEX option in superbuild.")
    list(REMOVE_ITEM CL_ARGS ${item})
  endif()
endforeach()

# Need to pass <INSTALL_DIR> to external project
# https://stackoverflow.com/questions/29723659/cmake-external-projects-command-seems-to-ignore-install-dir
list(APPEND CL_ARGS "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")

set(qmcpack_SOURCE_DIR "${CMAKE_SOURCE_DIR}/../")

include (ExternalProject)
ExternalProject_Add(
  qmcpack-real
  SOURCE_DIR ${qmcpack_SOURCE_DIR}
  CMAKE_ARGS "${CL_ARGS};-DQMC_COMPLEX=0"
)

ExternalProject_Add(
  qmcpack-complex
  SOURCE_DIR ${qmcpack_SOURCE_DIR}
  CMAKE_ARGS "${CL_ARGS};-DQMC_COMPLEX=1"
)

# Run tests for real and complex builds
# arguments passed through environment variable TESTS_ARGUMENTS
# https://stackoverflow.com/questions/28812533/how-to-pass-command-line-arguments-in-ctest-at-runtime
enable_testing()

ExternalProject_Get_Property(qmcpack-real binary_dir)
add_test(
    NAME qmcpack-real-test
    COMMAND ${CMAKE_COMMAND} -DTEST_EXECUTABLE=ctest -P ${qmcpack_SOURCE_DIR}/CMake/RunTestsSuperbuild.cmake
    WORKING_DIRECTORY ${binary_dir}
)

ExternalProject_Get_Property(qmcpack-complex binary_dir)
add_test(
    NAME qmcpack-complex-test
    COMMAND ${CMAKE_COMMAND} -DTEST_EXECUTABLE=ctest -P ${qmcpack_SOURCE_DIR}/CMake/RunTestsSuperbuild.cmake
    WORKING_DIRECTORY ${binary_dir}
)

# Install qmcpack and tools from qmcpack-real
ExternalProject_Get_Property(qmcpack-real install_dir)
install(DIRECTORY ${install_dir}/bin/
	DESTINATION bin
)

# Install qmcpack_complex from qmcpack-complex
ExternalProject_Get_Property(qmcpack-complex install_dir)
install(FILES ${install_dir}/bin/qmcpack_complex
        DESTINATION bin
)

install(FILES ${install_dir}/bin/qmcpack.settings
        DESTINATION bin
	RENAME qmcpack_complex.settings
)
