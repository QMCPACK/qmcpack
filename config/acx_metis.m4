dnl metis library is available at http://www-users.cs.umn.edu/~karypis/metis/
dnl
AC_DEFUN([ACX_METIS], [
acx_metis_ok=no

AC_ARG_WITH(metis,
        AC_HELP_STRING([--with-metis=DIR], 
                       [use metis library in DIR.]),
[
  metis_dir=$withval
  metis_include="$withval/include"
  metis_libs="-L$withval/lib -lmetis"
],
[])
                                                                                              
AC_ARG_WITH(metis-libs,
  AC_HELP_STRING([--with-metis-libs=<lib>],
                 [Use metis-related libraries <lib>]),
  [
    metis_libs="$withval"
  ],
  [])

dnl if test x"$with_metis" != x; then
if test x"$metis_dir" != x; then
  acx_metis_ok=yes
  METIS_LIBS="$metis_libs"
  if test -f "$metis_include/metis.h" ; then
    CPPFLAGS="$CPPFLAGS -I$metis_include"
  else
    acx_metis_ok=no
  fi
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_metis_ok" = xyes; then
  AC_SUBST(METIS_LIBS)
  ifelse([$1],,AC_DEFINE(HAVE_LIBMETIS,1,[Define if you have METIS library.]),[$1])
        :
else
  acx_metis_ok=no
  $2
fi
])dnl ACX_METIS
