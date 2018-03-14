# this new link path policy has fixed ldpath problems
# experienced with out of dyld_load_path libraries on
# OSX without bad effect on other platforms.
cmake_policy(SET CMP0003 NEW)

# in this directory and below remove the -DNDEBUG flag
string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
message( "CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE}" )
add_subdirectory(src)
