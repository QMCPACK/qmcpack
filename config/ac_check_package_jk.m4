dnl @synopsis AC_caolan_CHECK_PACKAGE(PACKAGE, FUNCTION, LIBRARY , HEADERFILE [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Provides --with-PACKAGE, --with-PACKAGE-include and --with-PACKAGE-libdir
dnl options to configure. Supports the now standard --with-PACKAGE=DIR 
dnl approach where the package's include dir and lib dir are underneath DIR,
dnl but also allows the include and lib directories to be specified seperately
dnl
dnl adds the extra -Ipath to CFLAGS if needed 
dnl adds extra -Lpath to LD_FLAGS if needed
dnl searches for the FUNCTION in the LIBRARY with 
dnl AC_CHECK_LIBRARY and thus adds the lib to LIBS
dnl
dnl defines HAVE_PKG_PACKAGE if it is found, (where PACKAGE in the 
dnl HAVE_PKG_PACKAGE is replaced with the actual first parameter passed)
dnl note that autoheader will complain of not having the HAVE_PKG_PACKAGE and you 
dnl will have to add it to acconfig.h manually
dnl
dnl @version $Id$
dnl @author Caolan McNamara <caolan@skynet.ie>
dnl
dnl with fixes from...
dnl Alexandre Duret-Lutz <duret_g@lrde.epita.fr>
dnl
dnl Does not work properly with sprng or hdf5
dnl with fixed from Jeongnim Kim 
dnl --with-pkg=dir dir/include/pkg is used if found

AC_DEFUN([AC_jk_CHECK_PACKAGE],
[

AC_ARG_WITH($1,
[  --with-$1[=DIR]	root directory of $1 installation],
with_$1=$withval 
if test x"${with_$1}" != x; then
   if test -d "$withval/include/$1" ; then
	$1_include="$withval/include/$1" 
   else
	$1_include="$withval/include" 
   fi
	$1_libdir="$withval/lib"
fi
)

AC_ARG_WITH($1-include,
[  --with-$1-include=DIR        specify exact include dir for $1 headers],
$1_include="$withval")

AC_ARG_WITH($1-libdir,
[  --with-$1-libdir=DIR        specify exact library dir for $1 library
  --without-$1        disables $1 usage completely], 
$1_libdir="$withval")


if test x"${$1_libdir}" != x || test x"${$1_include}" != x; then
  OLD_LIBS=$LIBS
  OLD_LDFLAGS=$LDFLAGS
  OLD_CFLAGS=$CFLAGS
  OLD_CPPFLAGS=$CPPFLAGS
  OLD_CXXFLAGS=$CXXFLAGS
dnl backup flags to go back if fails to find the library
  no_good=no

  if test "${$1_libdir}" ; then
	LDFLAGS="$LDFLAGS -L${$1_libdir}"
  fi
  if test "${$1_include}" ; then
        CPPFLAGS="$CPPFLAGS -I${$1_include}"
	CFLAGS="$CFLAGS -I${$1_include}"
	CXXFLAGS="$CXXFLAGS -I${$1_include}"
  fi
dnl will set HAVE_LIBUP[$1]
   AC_CHECK_LIB($3,$2,,no_good=yes)
   AC_CHECK_HEADER($4,,no_good=yes)

  if test "$no_good" = yes; then
dnl	broken
	ifelse([$6], , , [$6])
	LIBS=$OLD_LIBS
	LDFLAGS=$OLD_LDFLAGS
	CPPFLAGS=$OLD_CPPFLAGS
	CFLAGS=$OLD_CFLAGS
	CXXFLAGS=$OLD_CXXFLAGS
  else
dnl	fixed
	ifelse([$5], , , [$5])
	AC_DEFINE(HAVE_PKG_$1)
  fi
fi
])
