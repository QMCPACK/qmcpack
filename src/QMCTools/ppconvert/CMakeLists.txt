# In this directory subtree, remove the -DNDEBUG flag from build configs that add it
# C-style assert was used to check function success. For example,
#   assert (parser.ReadWord(psp_name));
# Hence all the assertion must be turned on
if(CMAKE_CXX_FLAGS_RELWITHDEBINFO)
  string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
endif()
if(CMAKE_CXX_FLAGS_RELEASE)
  string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
endif()

string(REPLACE "-ffast-math" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})

if(CMAKE_CROSSCOMPILING)
  message(STATUS "Cannot run std::isnan tests when cross compiling. Assume std::isnan not working. "
                 "Use -DQMC_ISNAN=ON to override.")
  set(ISNAN_WORKS FALSE)
else()
  # check if std::isnan works or not
  # this check may fail when cross-compiling even if std::isnan actually functions properly.
  include(CheckCXXSourceRuns)
  check_cxx_source_runs(
"
#include <cmath>
int main()
{
  if (std::isnan(std::sqrt(-1.)))
    return 0;
  return 1;
}
"
  ISNAN_WORKS)
endif()

option(QMC_ISNAN "Enable Code paths relying std::isnan" ${ISNAN_WORKS})

# ppconvert only works when std::isnan works.
if(QMC_ISNAN)
  message(STATUS "std::isnan dependent code paths enabled. ppconvert enabled.")
  add_subdirectory(src)
  add_subdirectory(test)
else()
  message(STATUS "std::isnan dependent code paths disabled. ppconvert disabled.")
endif()
