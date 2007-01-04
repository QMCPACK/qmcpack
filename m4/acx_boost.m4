dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_boost.html
dnl
AC_DEFUN([ACX_BOOST], [
acx_boost_ok=no

AC_ARG_WITH(boost,
        AC_HELP_STRING([--with-boost=DIR], 
                       [use boost library in DIR.]),
[
  boost_dir=$withval
  boost_include="$withval"
],
[])
                                                                                              
AC_ARG_WITH(boost-libs,
  AC_HELP_STRING([--with-boost-libs=<lib>]),
  [
    boost_libs="$withval"
  ],
  [])

if test x"$boost_dir" = x; then
  boost_dir_env="${BOOST_HOME}"
  if test x"$boost_dir_env" != x; then
    boost_dir=$boost_dir_env
    boost_include="$boost_dir_env"
    boost_libs="$boost_dir_env"
  fi
fi

dnl if test x"$with_boost" != x; then
if test x"$boost_dir" != x; then
  acx_boost_ok=yes
  if test -f "$boost_include/boost/config.hpp" ; then
    CPPFLAGS="$CPPFLAGS -I$boost_include"
  else
    acx_boost_ok=no
  fi
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_boost_ok" = xyes; then
  ifelse([$1],,AC_DEFINE(HAVE_LIBBOOST,1,[Define if you have BOOST library.]),[$1])
        :
else
  acx_boost_ok=no
  $2
fi
])dnl ACX_BOOST
