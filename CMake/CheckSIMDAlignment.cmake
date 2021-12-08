# Check if AVX512 is activated in the compilation
# Since cross-compiling is not unusual on HPC systems (Cray),
# try_compile is robust against
try_compile(CXX_COMPILER_HAVE_AVX512_MACRO ${CMAKE_BINARY_DIR} ${PROJECT_CMAKE}/try_compile_sources/check_AVX512.cpp
            CMAKE_FLAGS "${CMAKE_CXX_FLAGS}")

if(CXX_COMPILER_HAVE_AVX512_MACRO)
  set(default_alignment 64)
else()
  set(default_alignment 32)
endif()
