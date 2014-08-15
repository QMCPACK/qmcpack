#ADIOS_INCLUDES
#ADIOS_LIBRARIES
#ADIOS_FOUND
#ADIOS_DEFINITIONS

SET(ADIOS_HOME $ENV{ADIOS_HOME})
SET(ADIOS_CONFIG adios_config)
MESSAGE(STATUS "Try ADIOS_HOME="${ADIOS_HOME})

#Try to find adios config
IF (ADIOS_HOME)
  MESSAGE(STATUS "Found ADIOS_HOME="${ADIOS_HOME})
  FIND_PATH(ADIOS_INCLUDES adios.h ${ADIOS_HOME}/include)
  SET(ADIOS_FOUND 1)
  SET(ADIOS_CONFIG ${ADIOS_HOME}/bin/adios_config)   
  SET(ADIOS_LIB_PATHS ${ADIOS_HOME}/lib)
ELSE(ADIOS_HOME)
  #check for adios_config
  EXECUTE_PROCESS(COMMAND ${ADIOS_CONFIG} -d
                  OUTPUT_VARIABLE ADIOS_HOME
                  RESULT_VARIABLE ADIOS_FOUND_CONFIG)
  #Use adios config to find the correct flags
  IF (ADIOS_FOUND_CONFIG EQUAL 0)
    FIND_PATH(ADIOS_INCLUDES adios.h ${ADIOS_HOME}/include)
    SET(ADIOS_FOUND 1)
  ENDIF(ADIOS_FOUND_CONFIG EQUAL 0)
ENDIF(ADIOS_HOME)

SET(ADIOS_LIBRARIES "")


IF(ADIOS_FOUND)
  SET(FOUND_LIBRARYS 1)
  #Configure ADIOS with the correct flags
  EXECUTE_PROCESS(COMMAND ${ADIOS_CONFIG} -l  
    COMMAND sed s/-l//g  
    COMMAND sed s/-L[^\ ]*//g
    COMMAND sed s/-I[^\ ]*//g
    OUTPUT_VARIABLE ADIOS_CON_OUT)
 
  STRING(REPLACE "\n" "" ADIOS_CON_OUT "${ADIOS_CON_OUT}")
  STRING(REPLACE " " ";" ADIOS_CON_OUT "${ADIOS_CON_OUT}")
  #MESSAGE(STATUS ${ADIOS_CON_OUT})
  
  FOREACH (lib ${ADIOS_CON_OUT})
    #convert the library name to uppercase
    STRING(TOUPPER ${lib} NAME)
    IF(${NAME} STREQUAL "PTHREAD")
    ELSE(${NAME} STREQUAL "PTHREAD")
    #check the path for the library
    FIND_LIBRARY(${NAME}_LIBRARY ${lib} ${ADIOS_LIB_PATHS})
    #If we found the library then add it adios_libraries
    IF (${NAME}_LIBRARY)
      SET(ADIOS_LIBRARIES ${ADIOS_LIBRARIES} ${${NAME}_LIBRARY})
    #Call error if we can't find one of the libraries
    ELSE(${NAME}_LIBRARY)
      SET(FOUND_LIBRARYS 0)
      MESSAGE(STATUS "Try setting: " ${NAME}_LIBRARY)
    ENDIF(${NAME}_LIBRARY)
    ENDIF(${NAME} STREQUAL "PTHREAD")
  ENDFOREACH(lib)
  IF (${FOUND_LIBRARYS} EQUAL 0)
  	SET(ADIOS_FOUND 0)
    MESSAGE(FATAL_ERROR "ADIOS Dependencies were not all satisfied")
  ENDIF()
ENDIF(ADIOS_FOUND)


