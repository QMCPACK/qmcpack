#
# this module look for blitz++ (http://www.oonumerics.org/blitz) support
# it will define the following values
#
# BLITZ_INCLUDE_PATH = where blitz/blitz.h can be found
#

FIND_PATH(BLITZ_INCLUDE_PATH
    blitz/blitz.h
    /u/ac/esler/lib/blitz/include
    /home/common/lib/blitz-0.6-GCC
    /usr/apps/tools/blitz
    /usr/include
    /opt/include
    /usr/local/include
    ${AUXPACKAGES}/Utilities/blitz
)
INCLUDE_DIRECTORIES(${BLITZ_INCLUDE_PATH})

