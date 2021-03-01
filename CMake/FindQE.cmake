# Locate QE via the patched pw2qmcpack.x
# Require both to be present to consider QE_FOUND
# Take QE_BIN as hint for location

FIND_PATH(QE_PW_DIR pw.x HINTS ${QE_BIN})
FIND_PATH(QE_PW2Q_DIR pw2qmcpack.x HINTS ${QE_BIN})

SET(QE_FOUND FALSE)
IF(QE_PW2Q_DIR AND QE_PW_DIR)
  IF ( NOT (QE_PW2Q_DIR STREQUAL QE_PW_DIR) )
    MESSAGE(WARNING "Found pw.x and pw2qmcpack.x in different locations, ${QE_PW_DIR} and ${QE_PW2Q_DIR}, verify this is intentional.")
  ENDIF()
  #MESSAGE(STATUS "QE_PW2Q_DIR=${QE_PW2Q_DIR}")
  #MESSAGE(STATUS "QE_PW_DIR=${QE_PW_DIR}")
  SET(QE_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
   QE_PW2Q_DIR
   QE_PW_DIR
   QE_FOUND
)
