##########################################################################
# Boost Version information
# 1.34.x supports cmake build.
##########################################################################
set(BOOST_VERSION_MAJOR 1)
set(BOOST_VERSION_MINOR 33)
set(BOOST_VERSION_SUBMINOR 0)
set(BOOST_VERSION "${BOOST_VERSION_MAJOR}.${BOOST_VERSION_MINOR}.${BOOST_VERSION_SUBMINOR}")

SUBDIRS(serialization)
