dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_hdf5.html
dnl
AC_DEFUN([ACX_HDF5], [
acx_hdf5_ok=no

AC_ARG_WITH(hdf5,
        AC_HELP_STRING([--with-hdf5=DIR], 
                       [use hdf5 library in DIR.]),
[
  hdf5_dir=$withval
  hdf5_include="$withval/include"
  hdf5_libs="-L$withval/lib -lhdf5"
],
[])
                                                                                              
AC_ARG_WITH(hdf5-libs,
  AC_HELP_STRING([--with-hdf5-libs=<lib>],
                 [Use hdf5-related libraries <lib>]),
  [
    hdf5_libs="$withval"
  ],
  [])

dnl if test x"$with_hdf5" != x; then
if test x"$hdf5_dir" != x; then
  acx_hdf5_ok=yes
  HDF5_LIBS="$hdf5_libs"
  if test -f "$hdf5_include/hdf5.h" ; then
    CPPFLAGS="$CPPFLAGS -I$hdf5_include"
  else
    acx_hdf5_ok=no
  fi
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_hdf5_ok" = xyes; then
  AC_SUBST(HDF5_LIBS)
  ifelse([$1],,AC_DEFINE(HAVE_LIBHDF5,1,[Define if you have HDF5 library.]),[$1])
        :
else
  acx_hdf5_ok=no
  $2
fi
])dnl ACX_HDF5
