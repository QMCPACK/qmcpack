dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_sprng.html
dnl
AC_DEFUN([ACX_SPRNG], [
acx_sprng_ok=no

AC_ARG_WITH(sprng,
        AC_HELP_STRING([--with-sprng=DIR], [use sprng library in DIR.]),
[
  sprng_dir=$withval
  sprng_include="$withval/include"
  sprng_libs="-L$withval/lib -lsprng"
],
[])

AC_ARG_WITH(sprng-libs,
  AC_HELP_STRING([--with-sprng-libs=<lib>],
                 [Use sprng-related libraries <lib>]),
  [
    sprng_libs="$withval"
  ],
  [])

if test x"$sprng_dir" != x; then
  acx_sprng_ok=yes
  SPRNG_LIBS="$sprng_libs"
  if test -f "$sprng_include/sprng.h" ; then
    CPPFLAGS="$CPPFLAGS -I$sprng_include"
  else
    acx_sprng_ok=no
  fi
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_sprng_ok" = xyes; then
  AC_SUBST(SPRNG_LIBS)
  ifelse([$1],,AC_DEFINE(HAVE_LIBSPRNG,1,[Define if you have SPRNG library.]),[$1])
        :
else
  acx_sprng_ok=no
  $2
fi
])dnl ACX_SPRNG
