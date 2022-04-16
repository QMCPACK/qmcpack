# Locate rmg-cpu
# Take RMG_BIN as hint for location

find_program(RMG_CPU_EXE rmg-cpu HINTS ${RMG_BIN})

set(RMG_FOUND FALSE)
if(RMG_CPU_EXE)
  message(STATUS "RMG_CPU_EXE=${RMG_CPU_EXE}")
  set(RMG_FOUND TRUE)
endif()

mark_as_advanced(RMG_CPU_EXE RMG_FOUND)
