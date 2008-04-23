#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.61)
AC_INIT(qmcPACK, 0.2, esler@uiuc.edu, qmcPACK)
AM_INIT_AUTOMAKE(qmpPACK, 0.2)
AC_CONFIG_SRCDIR([src/ParticleTags.h])
AC_CONFIG_HEADER([src/config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL
PKG_PROG_PKG_CONFIG

# Optimal compiler flags
AX_CC_MAXOPT
AX_CXX_MAXOPT

AC_CHECK_LIB([m], [sqrt])
AC_MSG_CHECKING([LIBS = $LIBS])
# Checks for header files.
AC_HEADER_STDC

# Checks for language features
AC_C_VARARRAYS
AC_CXX_NAMESPACES

AC_MSG_CHECKING([optimal compiler flags])
AC_MSG_RESULT([${CFLAGS}])

###############################
# Checks for needed libraries #
###############################
AM_PATH_XML2([2.5.0])
CPPFLAGS="$CPPFLAGS $XML_CPPFLAGS"
LIBS="$LIBS $XML_LIBS"

AX_BOOST_BASE([1.34])

# Check if we want static-only binaries
AC_ARG_ENABLE(all-static, AC_HELP_STRING([--enable-all-static], [build static binaries (default=no)]))
if test "$enable_all_static" = "yes"; then
  ALL_STATIC="-all-static"
else
  ALL_STATIC=""
fi
AC_SUBST(ALL_STATIC)

#########################
# MPI and OpenMP checks #
#########################
AC_ARG_ENABLE(mpi, AC_HELP_STRING([--enable-mpi], [use MPI for message passing (default=no)]))
AM_CONDITIONAL(NEEDTOBUILDOOMPI,[test "$enable_mpi" = "yes"])
if test "$enable_mpi" = "yes"; then
  AC_DEFINE([HAVE_MPI],   [], [Define if compiling with MPI])
  AC_DEFINE([HAVE_OOMPI], [], [Define if using OOMPI]       )
fi
AC_ARG_WITH(mpi-cflags, AC_HELP_STRING([--with-mpi-cflags=],[preprocessor flags for compiling with mpi]))
AC_ARG_WITH(mpi-lib,    AC_HELP_STRING([--with-mpi-cflags=],[flags for linking with mpi]))


AC_ARG_ENABLE(openmp, AC_HELP_STRING([--enable-openmp], [use OpenMP for SMP parallelism (default=no)]))
if test "$enable_openmp" = "yes"; then
  AC_DEFINE([ENABLE_OPENMP], [], [Define if using OpenMP for threading.])
fi

AC_ARG_ENABLE([real]   , AC_HELP_STRING([--enable-real],    [build binary for real wave functions (default=yes)]),
              [],[enable_real="yes"])
AM_CONDITIONAL(WANT_REAL,[test "$enable_real" = "yes"])

AC_ARG_ENABLE([complex], AC_HELP_STRING([--enable-complex], [build binary for complex wave functions (default=yes)]),
              [], [enable_complex="yes"])
AM_CONDITIONAL(WANT_COMPLEX,[test "$enable_complex" = "yes"])


###########################
# OHMMS-related arguments #
###########################
AC_ARG_ENABLE([single-precision], 
              AC_HELP_STRING([--enable-single-precision], [build with 32-bit floats (default=no)]))
if test "$enable_single_precision" = "yes"; then
  AC_DEFINE([OHMMS_PRECISION], [float], [default floating point type])
else
  AC_DEFINE([OHMMS_PRECISION], [double], [default float point type])
fi

AC_ARG_WITH([dim],AC_HELP_STRING([--with-dim=n],[dimensionality of space (default=3)]))
if test "x$with_dim" = "x"; then
  with_dim=3
fi
AC_DEFINE([OHMMS_DIM],[3],[Dimensionality of space])
AC_DEFINE([OHMMS_INDEXTYPE],[int],[Define integer for the index type.])


# Profiling and debugging
AC_ARG_ENABLE(profile, AC_HELP_STRING([--enable-profile],[instrument code with profiling information]))
if test "x$enable_profile" = "xyes"; then
   CFLAGS="-pg $CFLAGS"
   FFLAGS="-pg $FFLAGS"
   LDFLAGS="-pg $LDFLAGS"
fi

AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug],[enable code for debugging checks]))
if test "x$enable_debug" = "xyes"; then
   CFLAGS="-DDEBUG $CFLAGS"
   FFLAGS="-DDEBUG $FFLAGS"
fi

AC_F77_WRAPPERS
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_SUBST(PKGDATADEF)

AC_CONFIG_FILES([ Makefile                  \
		  src/Makefile              \
                  qmcpack.pc  ])
AC_OUTPUT
