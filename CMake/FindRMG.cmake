# Locate rmg-cpu 
# Take RMG_BIN as hint for location

find_path(RMG_CPU_DIR rmg-cpu HINTS ${RMG_BIN})

set(RMG_FOUND FALSE)
if(RMG_CPU_DIR)
  MESSAGE(STATUS "RMG_CPU_DIR=${RMG_CPU_DIR}")
  set(RMG_FOUND TRUE)
endif()

mark_as_advanced(RMG_CPU_DIR RMG_FOUND)
