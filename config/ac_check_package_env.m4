dnl @synopsis AC_CHECK_PACKAGE_ENV(PACKAGE,PACKAGE_ENV_NAME [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Provides --enable-PACKAGE
dnl options to configure using environment variables.
dnl
dnl adds the extra include to CFLAGS, CPPFLAGS and CXXFLAGS
dnl adds the extra library to LDFLAGS 
dnl searches PACKAGE, PACKAGE_HOME, PACKAGE_LIBS and PACLAGE_INCLUDE
dnl
dnl defines HAVE_PKG_PACKAGE if it is found, (where PACKAGE in the 
dnl HAVE_PKG_PACKAGE is replaced with the actual first parameter passed)
dnl note that autoheader will complain of not having the HAVE_PKG_PACKAGE and you 
dnl will have to add it to acconfig.h manually
dnl
dnl @version $Id$
dnl @author Jeongnim Kim <jnkim@ncsa.uiuc.edu>
dnl
dnl with fixes from...
dnl 
dnl

AC_DEFUN([AC_CHECK_PACKAGE_ENV],
[

AC_ARG_ENABLE($1,
  AC_HELP_STRING([--enable-$1],[search $1 installation with $2 env (default=no)]),
  [enable_$1=$withval], 
  [enable_$1=no])

  if test x"${enable_$1}" != x"no"; then
    echo "checking for $1 with environment variables $2, $2_LIBS and $2_INCLUDE"
    $1_root="${$2}"
    $1_lib_env="${$2_LIBS}"
    $1_inc_env="${$2_INCLUDE}"
    if test x"${$1_root}" != x; then
      if test x"${$1_inc_env}" = x; then
        $1_inc_env="-I${$1_root}/include"
      fi 
      if test x"${$1_lib_env}" = x; then
        $2_LIBS=" -L${$1_root}/lib -l$1"
      fi 
      ifelse([$3], , , [$3])
      CPPFLAGS="$CPPFLAGS ${$1_inc_env}"
      CFLAGS="$CFLAGS ${$1_inc_env}"
      CXXFLAGS="$CXXFLAGS ${$1_inc_env}"
      AC_SUBST(${$2_LIBS})
    else
      ifelse([$4], , , [$4])
    fi
#      AC_SUBST(${$2_LIBS})
#      AC_DEFINE_UNQUOTED(HAVE_LIB$2, 1, [Define to 1 if you have $1 library])
#      echo "*** HAVE_LIB$2 is enabled"
#      echo "*** $2_LIBS is set to ${$2_LIBS}"
#      echo "*** ${$1_inc_env} is added to CPPFLAGS, CFLAGS and CXXFLAGS"
#    fi
#  else
#      AC_DEFINE_UNQUOTED(HAVE_LIB$2, 0, [Undefine to 1 if you don't have $1 library])
  fi
])
