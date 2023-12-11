# Locate QE via the patched pw2qmcpack.x
# Require both to be present to consider QE_FOUND
# Take QE_BIN as hint for location

find_path(QE_PW_DIR pw.x HINTS ${QE_BIN})
find_path(QE_PW2Q_DIR pw2qmcpack.x HINTS ${QE_BIN})

set(QE_FOUND FALSE)
if(QE_PW2Q_DIR AND QE_PW_DIR)
  if(NOT (QE_PW2Q_DIR STREQUAL QE_PW_DIR))
    message(
      WARNING
        "Found pw.x and pw2qmcpack.x in different locations, ${QE_PW_DIR} and ${QE_PW2Q_DIR}, verify this is intentional."
    )
  endif()
  #MESSAGE(STATUS "QE_PW2Q_DIR=${QE_PW2Q_DIR}")
  #MESSAGE(STATUS "QE_PW_DIR=${QE_PW_DIR}")
  set(QE_FOUND TRUE)
endif()

mark_as_advanced(QE_PW2Q_DIR QE_PW_DIR QE_FOUND)
