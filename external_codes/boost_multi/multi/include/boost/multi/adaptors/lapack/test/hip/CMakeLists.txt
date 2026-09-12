cmake_minimum_required(VERSION 3.24)

project(boost-multi-adaptors-lapack-test-hip LANGUAGES CXX HIP)

# GPU (hipSOLVER) tests for the lapack adaptor; unlike the CPU tests, these need neither
# a CPU LAPACK nor ENABLE_CIRCLE, so they are wired directly from the root under ENABLE_HIP

find_package(Boost REQUIRED COMPONENTS headers)

find_package(hipsolver REQUIRED)
find_package(rocthrust REQUIRED)

enable_testing()
include(CTest)

add_executable(lapack_cusolver.hip.x cusolver.hip)
target_compile_definitions(lapack_cusolver.hip.x PRIVATE MULTI_USE_HIP)  # redundant with lapack/hipsolver.hpp and the CI flags, but keeps stray TUs consistent
target_link_libraries(lapack_cusolver.hip.x PRIVATE multi roc::hipsolver roc::rocthrust Boost::headers)

add_test(NAME lapack_cusolver.hip.x COMMAND $<TARGET_FILE:lapack_cusolver.hip.x>)
