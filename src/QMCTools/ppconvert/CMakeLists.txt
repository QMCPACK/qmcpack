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

if(CMAKE_CXX_FLAGS MATCHES " -ffast-math")
  string(REPLACE " -ffast-math" " -fno-fast-math" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
endif()

if(CMAKE_CROSSCOMPILING)
  message(STATUS "Cannot run std::isnan tests when cross compiling. Setting ISNAN_WORKS to FALSE.")
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

option(ENABLE_PPCONVERT "Enable ppconvert" ${ISNAN_WORKS})

# ppconvert only works when std::isnan works.
if(ENABLE_PPCONVERT)
  message(STATUS "ppconvert enabled.")
  add_subdirectory(src)
  add_subdirectory(test)
else()
  message(STATUS "ppconvert disabled.")
endif()
