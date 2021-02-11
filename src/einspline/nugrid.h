/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef NUGRID_H
#define NUGRID_H

#include <stdbool.h>


typedef enum { LINEAR, GENERAL, CENTER, LOG } grid_type;

// Nonuniform grid base structure
typedef struct
{
  // public data
  grid_type code;
  double start, end;
  double* restrict points;
  int num_points;
  int (*reverse_map)(void *grid, double x);
} NUgrid;

#ifdef __cplusplus
extern "C"
#endif


typedef struct
{
  // public data
  grid_type code;
  double start, end;
  double* restrict points;
  int num_points;
  int (*reverse_map)(void *grid, double x);

  // private data
  double a, aInv, b, bInv, center, even_half;
  int half_points, odd_one;
  bool odd;
} center_grid;


typedef struct
{
  // public data
  grid_type code;
  double start, end;
  double* restrict points;
  int num_points;
  int (*reverse_map)(void *grid, double x);

  // private data
  double a, ainv, startinv;
} log_grid;


#ifdef __cplusplus
extern "C" {
#endif

  NUgrid*
  create_center_grid (double start, double end, double ratio,
                      int num_points);

  NUgrid*
  create_log_grid (double start, double end, int num_points);

  NUgrid*
  create_general_grid (double *points, int num_points);

  void
  destroy_grid (NUgrid *grid);

#ifdef __cplusplus
}
#endif
#endif
