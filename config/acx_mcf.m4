dnl mcf 
dnl
dnl
AC_DEFUN([ACX_MCF], [
acx_mcf_ok=no

AC_ARG_WITH(mcf,
        AC_HELP_STRING([--with-mcf=DIR], [use mcf library in DIR.]),
[
  mcf_dir=$withval
  mcf_include="$withval/include"
  mcf_libs="-L$withval/lib -lStormRT-v3.2.1"
],
[])

AC_ARG_WITH(mcf-libs,
  AC_HELP_STRING([--with-mcf-libs=<lib>],
                 [Use mcf-related libraries <lib>]),
  [
    mcf_libs="$withval"
  ],
  [])

if test x"$mcf_dir" != x; then
  acx_mcf_ok=yes
  MCF_LIBS="$mcf_libs"
  if test -f "$mcf_include/mcf.h" ; then
    CPPFLAGS="$CPPFLAGS -I$mcf_include"
    CXXFLAGS="${CXXFLAGS} -DSTORM_USEF -DUSE_STD -DSTORM_x86_LINUX"
  else
    acx_mcf_ok=no
  fi
fi

AM_CONDITIONAL(ENABLE_RTMRA, test x"$acx_mcf_ok" = x"yes")
# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_mcf_ok" = xyes; then
  AC_SUBST(MCF_LIBS)
  ifelse([$1],,AC_DEFINE(HAVE_LIBMCF,1,[Define if you have MCF library.]),[$1])
        :
else
  acx_mcf_ok=no
  $2
fi
])dnl ACX_MCF
