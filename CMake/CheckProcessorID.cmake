#
#checking processor name to pick optimization flags for linux and mac os x
#
IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
  exec_program(grep  
    ARGS model /proc/cpuinfo
    OUTPUT_VARIABLE CPUINFO_MODEL
    RETURN_VALUE CPUINFO
    )

  IF(CPUINFO_MODEL MATCHES "E5")
    SET(CPU_IDENTITY "E5xxx")
    MESSAGE("-- CPU Model Name Quad-core Xeon")
  ENDIF(CPUINFO_MODEL MATCHES "E5")

  IF(CPUINFO_MODEL MATCHES "E4")
    SET(CPU_IDENTITY "E4xxx")
    MESSAGE("-- CPU Model Name Core 2 Duo")
  ENDIF(CPUINFO_MODEL MATCHES "E4")

  IF(CPUINFO_MODEL MATCHES "Itanium")
    SET(CPU_IDENTITY "IA64")
    MESSAGE("-- CPU Model Name Itanium family")
  ENDIF(CPUINFO_MODEL MATCHES "Itanium")

ENDIF(CMAKE_SYSTEM_NAME MATCHES "Linux")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  exec_program(system_profiler
    ARGS -detailLevel basic
    OUTPUT_VARIABLE MACSYSINFO_MODEL
    RETURN_VALUE CPUINFO                                            
    )
  IF(MACSYSINFO_MODEL MATCHES "Core 2 Duo")
    SET(CPU_IDENTITY "E4xxx")
    MESSAGE("-- CPU Model Name Core 2 Duo ")
  ENDIF(MACSYSINFO_MODEL MATCHES "Core 2 Duo")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
