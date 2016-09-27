//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#include "multi_bspline.h"
#include "bspline.h"
#include "multi_nubspline.h"
#include "nubspline.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
  #include <omp.h>
#endif _OPENMP

double drand48();

inline double get_time()
{
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}

void
time_3d_real_double_omp()
{
  // int avail = numa_available();
#ifdef _OPENMP
  int nthr = omp_get_max_threads();
#else
  int nthr = 1;
#endif
  // int nnodes = numa_max_node();
  // fprintf (stderr, "Performing test with %d NUMA nodes.\n",
  // 	   avail, nnodes);
  // if (!nnodes)
  //   nnodes++;

  int nnodes = nthr;
  fprintf (stderr, "Using %d threads.\n", nnodes);

  int Nx=63; int Ny=61; int Nz = 69;
  int num_splines = 256;

  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 3.1; x_grid.end =  9.1; x_grid.num = Nx;
  y_grid.start = 8.7; y_grid.end = 12.7; y_grid.num = Ny;
  z_grid.start = 4.5; z_grid.end =  9.3; z_grid.num = Nz;

  BCtype_d xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = PERIODIC;
  yBC.lCode = yBC.rCode = PERIODIC;
  zBC.lCode = zBC.rCode = PERIODIC;

  // First, create splines the normal way
  UBspline_3d_d* norm_splines[num_splines];
  multi_UBspline_3d_d *multi_spline[nnodes];
  
  // First, create multispline
#pragma omp parallel for
  for (int node=0; node<nnodes; node++) 
  {
    // nodemask_t mask;
    // nodemask_zero(&mask);
    // nodemask_set (&mask, node);
    // numa_set_membind (&mask);
    multi_spline[node] = create_multi_UBspline_3d_d 
      (x_grid, y_grid, z_grid, xBC, yBC, zBC, num_splines);
  }

  double data[Nx*Ny*Nz];
  // Now, create normal splines and set multispline data
  for (int i=0; i<num_splines; i++) {
    for (int j=0; j<Nx*Ny*Nz; j++)
      data[j] = (drand48()-0.5);
    norm_splines[i] = create_UBspline_3d_d 
      (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
#pragma omp parallel for    
    for (int node=0; node<nnodes; node++) {
      // nodemask_t mask;
      // nodemask_zero(&mask);
      // nodemask_set (&mask, node);
      // numa_set_membind (&mask);
      set_multi_UBspline_3d_d (multi_spline[node], i, data);
    }
  }
  
  // Now, test random values
  double rand_start, rand_end, norm_start[nthr], norm_end[nthr], multi_start[nthr], multi_end[nthr];
  int num_vals = 10000;
  double multi_vals[nthr][num_splines], norm_vals[nthr][num_splines];
  double multi_grads[nthr][3*num_splines], norm_grads[nthr][3*num_splines];
  double multi_lapl[nthr][num_splines], norm_lapl[nthr][num_splines];
  double multi_hess[nthr][9*num_splines], norm_hess[nthr][9*num_splines];

  rand_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
  }
  rand_end = get_time();

  ///////////////////////
  // Check value routine  //
  ///////////////////////
  double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
  double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
  double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;

  int thr_per_node = nthr/nnodes;

#pragma omp parallel for
  for (int thr=0; thr<nthr; thr++) {
    int node = thr/thr_per_node;
    multi_start[thr] = get_time();
    for (int i=0; i<num_vals; i++) {
      double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end; 
      double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end; 
      double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end; 
      eval_multi_UBspline_3d_d (multi_spline[node], x, y, z, multi_vals[thr]);
    }
    multi_end[thr] = get_time();
  }

// #pragma omp parallel for
//   for (int thr=0; thr<nthr; thr++) {
//     norm_start[thr] = get_time();
//     for (int i=0; i<num_vals; i++) {
//       double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
//       double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
//       double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
//       for (int j=0; j<num_splines; j++)
// 	eval_UBspline_3d_d (norm_splines[j], x, y, z, &(norm_vals[thr][j]));
//     }
//     norm_end[thr] = get_time();
//   }
  
  double norm_avg=0.0, multi_avg=0.0;

  for (int thr=0; thr<nthr; thr++) {
    double norm_time   = (double)(norm_end[thr] - norm_start[thr] + rand_start - rand_end);
    double multi_time  = (double)(multi_end[thr] - multi_start[thr] + rand_start - rand_end);
    norm_avg += norm_time;
    multi_avg += multi_time;
  }
  norm_avg  /= nthr;
  multi_avg /= nthr;
  double norm_speed  = (double) num_vals*num_splines / norm_avg;
  double multi_speed = (double) num_vals*num_splines / multi_avg;

  // fprintf (stderr, "Normal value speed = %13.3f evaluations per second.\n", 
  // 	   norm_speed);
  fprintf (stderr, "Multi  value speed = %13.3f evaluations per second.\n", 
  	   multi_speed);
  fprintf (stderr, "Aggregate bandwidth = %1.3f GB/s per socket\n", multi_speed * 64.0*8.0 * 8 * 1.0e-9);

  
  ///////////////////////
  // Check VGH routine //
  ///////////////////////
  #pragma omp parallel for
  for (int thr=0; thr<nthr; thr++) {
    int node = thr/thr_per_node;
    multi_start[thr] = get_time();
    for (int i=0; i<num_vals; i++) {
      double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
      double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
      double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
      eval_multi_UBspline_3d_d_vgh 
	(multi_spline[node], x, y, z,  multi_vals[thr], 
	 multi_grads[thr], multi_hess[thr]);
    }
    multi_end[thr] = get_time();
  }

// #pragma omp parallel for
//   for (int thr=0; thr<nthr; thr++) {
//     norm_start[thr] = get_time();
//     for (int i=0; i<num_vals; i++) {
//       double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
//       double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
//       double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
//       for (int j=0; j<num_splines; j++)
// 	eval_UBspline_3d_d_vgh (norm_splines[j], x, y, z, &(norm_vals[thr][j]),
// 				&(norm_grads[thr][3*j]), &(norm_hess[thr][9*j]));
//     }
//     norm_end[thr] = get_time();
//   }

  norm_avg = multi_avg = 0.0;
  for (int thr=0; thr<nthr; thr++) {
    double norm_time   = (double)(norm_end[thr] - norm_start[thr] + rand_start - rand_end);
    double multi_time  = (double)(multi_end[thr] - multi_start[thr] + rand_start - rand_end);
    norm_avg += norm_time;
    multi_avg += multi_time;
  }
  norm_avg  /= nthr;
  multi_avg /= nthr;
  norm_speed  = (double) num_vals*num_splines / norm_avg;
  multi_speed = (double) num_vals*num_splines / multi_avg;

//   fprintf (stderr, "Normal VGH   speed = %13.3f evaluations per second.\n", 
// 	   norm_speed);
  fprintf (stderr, "Multi  VGH   speed = %13.3f evaluations per second.\n", 
	   multi_speed);
  fprintf (stderr, "%1.3f GFLOPS per socket\n", multi_speed * 64.0*2.0*10.0 * 8 * 1.0e-9);



//   destroy_Bspline (multi_spline);
//   for (int i=0; i<num_splines; i++)
//     destroy_Bspline(norm_splines[i]); 
}


void
time_3d_complex_double_omp()
{
#ifdef _OPENMP
  int nthr = omp_get_max_threads();
#else
  int nthr = 1;
#endif
  int nnodes = nthr;
  fprintf (stderr, "Using %d threads.\n", nthr);

  int Nx=32; int Ny=32; int Nz = 32;
  int num_splines = 256;

  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 3.1; x_grid.end =  9.1; x_grid.num = Nx;
  y_grid.start = 8.7; y_grid.end = 12.7; y_grid.num = Ny;
  z_grid.start = 4.5; z_grid.end =  9.3; z_grid.num = Nz;

  BCtype_z xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = PERIODIC;
  yBC.lCode = yBC.rCode = PERIODIC;
  zBC.lCode = zBC.rCode = PERIODIC;

  // First, create splines the normal way
  UBspline_3d_z* norm_splines[num_splines];
  multi_UBspline_3d_z *multi_spline[nthr];
  
  // First, create multispline
#pragma omp parallel for
  for (int node=0; node<nthr; node++) 
  {
    // nodemask_t mask;
    // nodemask_zero(&mask);
    // nodemask_set (&mask, node);
    // numa_set_membind (&mask);
    multi_spline[node] = create_multi_UBspline_3d_z
      (x_grid, y_grid, z_grid, xBC, yBC, zBC, num_splines);
  }

  double data[Nx*Ny*Nz*2];
  // Now, create normal splines and set multispline data
  for (int i=0; i<num_splines; i++) {
    for (int j=0; j<Nx*Ny*Nz; j++)
      data[j] = (drand48()-0.5);
    norm_splines[i] = create_UBspline_3d_z
      (x_grid, y_grid, z_grid, xBC, yBC, zBC, (complex_double*)data);
#pragma omp parallel for    
    for (int node=0; node<nthr; node++) {
      // nodemask_t mask;
      // nodemask_zero(&mask);
      // nodemask_set (&mask, node);
      // numa_set_membind (&mask);
      set_multi_UBspline_3d_z (multi_spline[node], i, data);
    }
  }
  
  // Now, test random values
  double rand_start, rand_end, norm_start[nthr], norm_end[nthr], multi_start[nthr], multi_end[nthr];
  int num_vals = 10000;
  complex_double multi_vals[nthr][num_splines], norm_vals[nthr][num_splines];
  complex_double multi_grads[nthr][3*num_splines], norm_grads[nthr][3*num_splines];
  complex_double multi_lapl[nthr][num_splines], norm_lapl[nthr][num_splines];
  complex_double multi_hess[nthr][9*num_splines], norm_hess[nthr][9*num_splines];

  rand_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
  }
  rand_end = get_time();

  ///////////////////////
  // Check value routine  //
  ///////////////////////
  double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
  double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
  double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;

  int thr_per_node = nthr/nthr;

#pragma omp parallel for
  for (int thr=0; thr<nthr; thr++) {
    int node = thr/thr_per_node;
    multi_start[thr] = get_time();
    for (int i=0; i<num_vals; i++) {
      double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end; 
      double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end; 
      double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end; 
      eval_multi_UBspline_3d_z (multi_spline[node], x, y, z, multi_vals[thr]);
    }
    multi_end[thr] = get_time();
  }

// #pragma omp parallel for
//   for (int thr=0; thr<nthr; thr++) {
//     norm_start[thr] = get_time();
//     for (int i=0; i<num_vals; i++) {
//       double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
//       double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
//       double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
//       for (int j=0; j<num_splines; j++)
// 	eval_UBspline_3d_z (norm_splines[j], x, y, z, &(norm_vals[thr][j]));
//     }
//     norm_end[thr] = get_time();
//   }
  
  double norm_avg=0.0, multi_avg=0.0;

  for (int thr=0; thr<nthr; thr++) {
    double norm_time   = (double)(norm_end[thr] - norm_start[thr] + rand_start - rand_end);
    double multi_time  = (double)(multi_end[thr] - multi_start[thr] + rand_start - rand_end);
    norm_avg += norm_time;
    multi_avg += multi_time;
  }
  norm_avg  /= nthr;
  multi_avg /= nthr;
  double norm_speed  = (double) num_vals*num_splines / norm_avg;
  double multi_speed = (double) num_vals*num_splines / multi_avg;

  // fprintf (stderr, "Normal value speed = %13.3f evaluations per second.\n", 
  // 	   norm_speed);
  fprintf (stderr, "Multi  value speed = %13.3f evaluations per second.\n", 
  	   multi_speed);
  fprintf (stderr, "Aggregate bandwidth = %1.3f GB/s per socket\n", multi_speed * 64.0*16.0 * 8 * 1.0e-9);
  fprintf (stderr, "%1.3f GFLOPS per socket\n", multi_speed * 64.0*4.0 * 8 * 1.0e-9);

  
  ///////////////////////
  // Check VGH routine //
  ///////////////////////
#pragma omp parallel for
  for (int thr=0; thr<nthr; thr++) {
    int node = thr/thr_per_node;
    multi_start[thr] = get_time();
    for (int i=0; i<num_vals; i++) {
      double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
      double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
      double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
      eval_multi_UBspline_3d_z_vgh 
	(multi_spline[node], x, y, z,  multi_vals[thr], 
	 multi_grads[thr], multi_hess[thr]);
    }
    multi_end[thr] = get_time();
  }

// #pragma omp parallel for
//   for (int thr=0; thr<nthr; thr++) {
//     norm_start[thr] = get_time();
//     for (int i=0; i<num_vals; i++) {
//       double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
//       double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
//       double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
//       for (int j=0; j<num_splines; j++)
// 	eval_UBspline_3d_z_vgh (norm_splines[j], x, y, z, &(norm_vals[thr][j]),
// 				&(norm_grads[thr][3*j]), &(norm_hess[thr][9*j]));
//     }
//     norm_end[thr] = get_time();
//   }

  norm_avg = multi_avg = 0.0;
  for (int thr=0; thr<nthr; thr++) {
    double norm_time   = (double)(norm_end[thr] - norm_start[thr] + rand_start - rand_end);
    double multi_time  = (double)(multi_end[thr] - multi_start[thr] + rand_start - rand_end);
    norm_avg += norm_time;
    multi_avg += multi_time;
  }
  norm_avg  /= nthr;
  multi_avg /= nthr;
  norm_speed  = (double) num_vals*num_splines / norm_avg;
  multi_speed = (double) num_vals*num_splines / multi_avg;

//   fprintf (stderr, "Normal VGH   speed = %13.3f evaluations per second.\n", 
// 	   norm_speed);
  fprintf (stderr, "Multi  VGH   speed = %13.3f evaluations per second.\n", 
	   multi_speed);
  fprintf (stderr, "%1.3f GFLOPS per socket\n", multi_speed * 64.0*4.0*10.0 * 8 * 1.0e-9);


//   destroy_Bspline (multi_spline);
//   for (int i=0; i<num_splines; i++)
//     destroy_Bspline(norm_splines[i]); 
}


main()
{
  // fprintf (stderr, "Real:\n");
  // time_3d_real_double_omp();
  fprintf (stderr, "\nComplex:\n");
  time_3d_complex_double_omp();
}
