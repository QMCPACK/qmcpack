# in this directory and below remove the -DNDEBUG flag from build configs that add it
string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
string(REPLACE "-ffast-math" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
add_subdirectory(src)
add_subdirectory(test)
