dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_blitz.html
dnl
AC_DEFUN([ACX_BLITZ], [
acx_blitz_ok=no

AC_ARG_WITH(blitz,
        AC_HELP_STRING([--with-blitz=DIR], 
                       [use blitz library in DIR.]),
[
  blitz_dir=$withval
  blitz_include="$withval"
],
[])
                                                                                              
dnl if test x"$with_blitz" != x; then
if test x"$blitz_dir" != x; then
  acx_blitz_ok=yes
  if test -f "$blitz_include/blitz/blitz.h" ; then
    CPPFLAGS="$CPPFLAGS -I$blitz_include"
  else
    acx_blitz_ok=no
  fi
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blitz_ok" = xyes; then
  ifelse([$1],,AC_DEFINE(HAVE_LIBBLITZ,1,[Define if you have BLITZ library.]),[$1])
        :
else
  acx_blitz_ok=no
  $2
fi
])dnl ACX_BLITZ
