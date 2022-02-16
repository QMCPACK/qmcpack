# Einspline

Some of the files included in this QMCPACK einspline directory are based on the einspline library
released by Ken Esler http://einspline.sourceforge.net/ . einspline is a C library for the creation
and evaluation of interpolating cubic basis splines (B-splines) in 1, 2, and 3 dimensions.

On 2018-02-20 the original einspline library was released under the BSD-3-clause license
( https://opensource.org/licenses/BSD-3-Clause ), as noted at
https://sourceforge.net/p/einspline/code/443/ .


# Folder structure
9 test related files not coverged by unit tests.

The CPU part. Then following head files should be included by QMC subroutines on demand
```
bspline.h		single unifrom bspline
multi_bspline.h		multiple unifrom bspline
```

suffix
```
_structs.h		data structure
_create.h		initialization function declariation
_eval_s/d/c/z.h		evaluation function declariation
_eval_?_std/sse.cpp	the implementation selected by CMakeLists.txt
```
