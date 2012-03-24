/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//                                                                         //
//  This program is free software; you can redistribute it and/or modify   //
//  it under the terms of the GNU General Public License as published by   //
//  the Free Software Foundation; either version 2 of the License, or      //
//  (at your option) any later version.                                    //
//                                                                         //
//  This program is distributed in the hope that it will be useful,        //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of         //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          //
//  GNU General Public License for more details.                           //
//                                                                         //
//  You should have received a copy of the GNU General Public License      //
//  along with this program; if not, write to the Free Software            //
//  Foundation, Inc., 51 Franklin Street, Fifth Floor,                     //
//  Boston, MA  02110-1301  USA                                            //
/////////////////////////////////////////////////////////////////////////////

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
#endif

double drand48();

inline double get_time()
{
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}

inline double diff (double a, double b, double tol)
{
  if (fabs(a-b) > tol) 
    return 1;
  else
    return 0;
}

inline int
zdiff (complex_double a, complex_double b, double tol)
{
  double rdiff = fabs(creal(a) - creal(b));
  double idiff = fabs(cimag(a) - cimag(b));
  if (rdiff > tol || idiff > tol)
    return 1;
  else
    return 0;
}


// int 
// test_3d_double_all()
// {
//   int Nx=73; int Ny=91; int Nz = 29;
//   int num_splines = 128;

//   Ugrid x_grid, y_grid, z_grid;
//   x_grid.start = 3.1; x_grid.end =  9.1; x_grid.num = Nx;
//   y_grid.start = 8.7; y_grid.end = 12.7; y_grid.num = Ny;
//   z_grid.start = 4.5; z_grid.end =  9.3; z_grid.num = Nz;

//   BCtype_d xBC, yBC, zBC;
//   xBC.lCode = xBC.rCode = PERIODIC;
//   yBC.lCode = yBC.rCode = PERIODIC;
//   zBC.lCode = zBC.rCode = PERIODIC;

//   // First, create splines the normal way
//   UBspline_3d_d* norm_splines[num_splines];
//   multi_UBspline_3d_d *multi_spline;
  
//   // First, create multispline
//   multi_spline = create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC,
// 					     num_splines);

//   double data[Nx*Ny*Nz];
//   // Now, create normal splines and set multispline data
//   for (int i=0; i<num_splines; i++) {
//     for (int j=0; j<Nx*Ny*Nz; j++)
//       data[j] = (drand48()-0.5);// + (drand48()-0.5)*1.0i;
//     norm_splines[i] = create_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
//     set_multi_UBspline_3d_d (multi_spline, i, data);
//   }

// //   fprintf (stderr, "norm coef  = %1.14e + %1.14ei\n",
// // 	   creal(norm_splines[19]->coefs[227]),
// // 	   cimag(norm_splines[19]->coefs[227]));
// //   fprintf (stderr, "multi coef = %1.14e + %1.14ei\n",
// // 	   creal(multi_spline->coefs[19+227*multi_spline->z_stride]),
// // 	   cimag(multi_spline->coefs[19+227*multi_spline->z_stride]));
  
//   // Now, test random values
//   int num_vals = 100;
//   double multi_vals[num_splines], norm_vals[num_splines];
//   double multi_grads[3*num_splines], norm_grads[3*num_splines];
//   double multi_lapl[num_splines], norm_lapl[num_splines];
//   double multi_hess[9*num_splines], norm_hess[9*num_splines];
//   for (int i=0; i<num_vals; i++) {
//     double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
//     double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
//     double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    

//     ///////////////////////
//     // Check VG routine  //
//     ///////////////////////
//     eval_multi_UBspline_3d_d_vg (multi_spline, x, y, z, 
// 				  multi_vals, multi_grads);
//     for (int j=0; j<num_splines; j++)
//       eval_UBspline_3d_d_vg (norm_splines[j], x, y, z, &(norm_vals[j]),
// 			  &(norm_grads[3*j]));
//     for (int j=0; j<num_splines; j++) {
//       // Check value
//       if (diff(norm_vals[j], multi_vals[j], 1.0e-12))
// 	return -1;
      
//       // Check gradients
//       for (int n=0; n<3; n++) 
// 	if (diff (norm_grads[3*j+n], multi_grads[3*j+n], 1.0e-12))
// 	  return -2;
//     }


//     ///////////////////////
//     // Check VGL routine //
//     ///////////////////////
//     eval_multi_UBspline_3d_d_vgl (multi_spline, x, y, z, 
// 				  multi_vals, multi_grads, multi_lapl);
//     for (int j=0; j<num_splines; j++)
//       eval_UBspline_3d_d_vgl (norm_splines[j], x, y, z, &(norm_vals[j]),
// 			  &(norm_grads[3*j]), &(norm_lapl[j]));
//     for (int j=0; j<num_splines; j++) {
//       // Check value
//       if (diff(norm_vals[j], multi_vals[j], 1.0e-12))
// 	return -3;

//       // Check gradients
//       for (int n=0; n<3; n++) 
// 	if (diff (norm_grads[3*j+n], multi_grads[3*j+n], 1.0e-10))
// 	  return -4;

//       // Check laplacian
//       if (diff (norm_lapl[j], multi_lapl[j], 1.0e-10)) 
// 	return -5;
//     }


//     ///////////////////////
//     // Check VGH routine //
//     ///////////////////////
//     eval_multi_UBspline_3d_d_vgh (multi_spline, x, y, z, 
// 				  multi_vals, multi_grads, multi_hess);
//     for (int j=0; j<num_splines; j++)
//       eval_UBspline_3d_d_vgh (norm_splines[j], x, y, z, &(norm_vals[j]),
// 			  &(norm_grads[3*j]), &(norm_hess[9*j]));
//     for (int j=0; j<num_splines; j++) {
//       // Check value
//       if (diff(norm_vals[j], multi_vals[j], 1.0e-12))
// 	return -6;

//       // Check gradients
//       for (int n=0; n<3; n++) 
// 	if (diff (norm_grads[3*j+n], multi_grads[3*j+n], 1.0e-12)) 
// 	  return -7;

//       // Check hessian
//       for (int n=0; n<9; n++) 
// 	if (diff (norm_hess[9*j+n], multi_hess[9*j+n], 1.0e-10)) 
// 	  return -8;
//     }
//   }
//   return 0;
// }



//////////////////////////////////////////
// Single-precision real test functions //
//////////////////////////////////////////


/* void
time_3d_complex_float_all()
{
  int Nx=23; int Ny=21; int Nz = 29;
  int num_splines = 128;

  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 3.1; x_grid.end =  9.1; x_grid.num = Nx;
  y_grid.start = 8.7; y_grid.end = 12.7; y_grid.num = Ny;
  z_grid.start = 4.5; z_grid.end =  9.3; z_grid.num = Nz;

  BCtype_c xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = PERIODIC;
  yBC.lCode = yBC.rCode = PERIODIC;
  zBC.lCode = zBC.rCode = PERIODIC;

  // First, create splines the normal way                                       
  UBspline_3d_c* norm_splines[num_splines];
  multi_UBspline_3d_c *multi_spline;

  // First, create multispline                                                  
  multi_spline = create_multi_UBspline_3d_c (x_grid, y_grid, z_grid, xBC, yBC, 
					     zBC, num_splines);

  complex_float data[Nx*Ny*Nz];
  // Now, create normal splines and set multispline data                        
  for (int i=0; i<num_splines; i++) {
    for (int j=0; j<Nx*Ny*Nz; j++)
      data[j] = (drand48()-0.5);// + (drand48()-0.5)*1.0i;
    norm_splines[i] = create_UBspline_3d_c (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
    set_multi_UBspline_3d_c (multi_spline, i, data);
  }

  // Now, test random values                                                    
  int num_vals = 100000;
  complex_float multi_vals[num_splines], norm_vals[num_splines];
  complex_float multi_grads[3*num_splines], norm_grads[3*num_splines];
  complex_float multi_lapl[num_splines], norm_lapl[num_splines];
  complex_float multi_hess[9*num_splines], norm_hess[9*num_splines];

  double rand_start, rand_end, norm_start, norm_end, multi_start, multi_end;

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
  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    eval_multi_UBspline_3d_c (multi_spline, x, y, z, multi_vals);
  }
  multi_end = get_time();

  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_c (norm_splines[j], x, y, z, &(norm_vals[j]));
  }
  norm_end = get_time();

  double norm_time   = (double)(norm_end - norm_start + rand_start - rand_end);

  double multi_time  = (double)(multi_end - multi_start + rand_start - rand_end) ;
  double norm_speed  = (double) num_vals*num_splines / norm_time;
  double multi_speed = (double) num_vals*num_splines / multi_time;
  fprintf (stderr, "Normal value speed = %13.3f evaluations per second.\n", norm_speed);
  fprintf (stderr, "Multi  value speed = %13.3f evaluations per second.\n", multi_speed);
  ///////////////////////                                                       
  // Check VGH routine //                                                       
  ///////////////////////                                                       
  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    eval_multi_UBspline_3d_c_vgh (multi_spline, x, y, z,
                                  multi_vals, multi_grads, multi_hess);
  }
  multi_end = get_time();

  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_c_vgh (norm_splines[j], x, y, z, &(norm_vals[j]),
			      &(norm_grads[3*j]), &(norm_hess[9*j]));
  }
  norm_end = get_time();

  norm_time   = (double)(norm_end - norm_start + rand_start - rand_end) ;
  multi_time  = (double)(multi_end - multi_start + rand_start - rand_end) ;
  norm_speed  = (double) num_vals*num_splines / norm_time;
  multi_speed = (double) num_vals*num_splines / multi_time;
  fprintf (stderr, "Normal VGH   speed = %13.3f evaluations per second.\n", norm_speed);
  fprintf (stderr, "Multi  VGH   speed = %13.3f evaluations per second.\n", multi_speed);

  destroy_Bspline (multi_spline);
  for (int i=0; i<num_splines; i++)
    destroy_Bspline(norm_splines[i]);
}*/



void
time_3d_real_float_all()
{
  int Nx=23; int Ny=21; int Nz = 29;
  int num_splines = 128;

  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 3.1; x_grid.end =  9.1; x_grid.num = Nx;
  y_grid.start = 8.7; y_grid.end = 12.7; y_grid.num = Ny;
  z_grid.start = 4.5; z_grid.end =  9.3; z_grid.num = Nz;

  BCtype_s xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = PERIODIC;
  yBC.lCode = yBC.rCode = PERIODIC;
  zBC.lCode = zBC.rCode = PERIODIC;

  // First, create splines the normal way
  UBspline_3d_s* norm_splines[num_splines];
  multi_UBspline_3d_s *multi_spline;
  
  // First, create multispline
  multi_spline = create_multi_UBspline_3d_s 
    (x_grid, y_grid, z_grid, xBC, yBC, zBC, num_splines);

  float data[Nx*Ny*Nz];
  // Now, create normal splines and set multispline data
  for (int i=0; i<num_splines; i++) {
    for (int j=0; j<Nx*Ny*Nz; j++)
      data[j] = (drand48()-0.5);
    norm_splines[i] = create_UBspline_3d_s 
      (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
    set_multi_UBspline_3d_s (multi_spline, i, data);
  }
  
  // Now, test random values
  int num_vals = 100000;
  float multi_vals[num_splines], norm_vals[num_splines];
  float multi_grads[3*num_splines], norm_grads[3*num_splines];
  float multi_lapl[num_splines], norm_lapl[num_splines];
  float multi_hess[9*num_splines], norm_hess[9*num_splines];

  double rand_start, rand_end, norm_start, norm_end, multi_start, multi_end;

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
  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    eval_multi_UBspline_3d_s (multi_spline, x, y, z, multi_vals);
  }
  multi_end = get_time();

  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_s (norm_splines[j], x, y, z, &(norm_vals[j]));
  }
  norm_end = get_time();
  
  double norm_time   = (double)(norm_end - norm_start + rand_start - rand_end) ;
  double multi_time  = (double)(multi_end - multi_start + rand_start - rand_end) ;
  double norm_speed  = (double) num_vals*num_splines / norm_time;
  double multi_speed = (double) num_vals*num_splines / multi_time;
  fprintf (stderr, "Normal value speed = %13.3f evaluations per second.\n", 
	   norm_speed);
  fprintf (stderr, "Multi  value speed = %13.3f evaluations per second.\n", 
	   multi_speed);
  

  ///////////////////////
  // Check VGH routine //
  ///////////////////////
  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    eval_multi_UBspline_3d_s_vgh (multi_spline, x, y, z, 
				  multi_vals, multi_grads, multi_hess);
  }
  multi_end = get_time();

  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_s_vgh (norm_splines[j], x, y, z, &(norm_vals[j]),
			  &(norm_grads[3*j]), &(norm_hess[9*j]));
  }
  norm_end = get_time();

  norm_time   = (double)(norm_end - norm_start + rand_start - rand_end) ;
  multi_time  = (double)(multi_end - multi_start + rand_start - rand_end) ;
  norm_speed  = (double) num_vals*num_splines / norm_time;
  multi_speed = (double) num_vals*num_splines / multi_time;
  fprintf (stderr, "Normal VGH   speed = %13.3f evaluations per second.\n", 
	   norm_speed);
  fprintf (stderr, "Multi  VGH   speed = %13.3f evaluations per second.\n", 
	   multi_speed);

  destroy_Bspline (multi_spline);
  for (int i=0; i<num_splines; i++)
    destroy_Bspline(norm_splines[i]); 
}



#ifdef _OPENMP

#include <omp.h>
//#include <numa.h>

void
time_3d_real_double_omp()
{
  // int avail = numa_available();
  int nthr = omp_get_max_threads();
  // int nnodes = numa_max_node();
  // fprintf (stderr, "Performing test with %d NUMA nodes.\n",
  // 	   avail, nnodes);
  // if (!nnodes)
  //   nnodes++;

  int nnodes = omp_get_num_threads();

  int Nx=63; int Ny=61; int Nz = 69;
  int num_splines = 128;

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
  //#pragma omp parallel for
  for (int node=0; node<nnodes; node++) {
    // nodemask_t mask;
    // nodemask_zero(&mask);
    // nodemask_set (&mask, node);
    // numa_set_membind (&mask);
    // multi_spline[node] = create_multi_UBspline_3d_d 
    //   (x_grid, y_grid, z_grid, xBC, yBC, zBC, num_splines);
  }

//   double data[Nx*Ny*Nz];
//   // Now, create normal splines and set multispline data
//   for (int i=0; i<num_splines; i++) {
//     for (int j=0; j<Nx*Ny*Nz; j++)
//       data[j] = (drand48()-0.5);
//     norm_splines[i] = create_UBspline_3d_d 
//       (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
//     for (int node=0; node<nnodes; node++) {
//       nodemask_t mask;
//       nodemask_zero(&mask);
//       nodemask_set (&mask, node);
//       numa_set_membind (&mask);
//       set_multi_UBspline_3d_d (multi_spline[node], i, data);
//     }
//   }
  
//   // Now, test random values
//   double rand_start, rand_end, norm_start[nthr], norm_end[nthr], multi_start[nthr], multi_end[nthr];
//   int num_vals = 100000;
//   double multi_vals[nthr][num_splines], norm_vals[nthr][num_splines];
//   double multi_grads[nthr][3*num_splines], norm_grads[nthr][3*num_splines];
//   double multi_lapl[nthr][num_splines], norm_lapl[nthr][num_splines];
//   double multi_hess[nthr][9*num_splines], norm_hess[nthr][9*num_splines];

//   rand_start = omp_get_wtime();
//   for (int i=0; i<num_vals; i++) {
//     double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
//     double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
//     double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
//   }
//   rand_end = omp_get_wtime();

//   ///////////////////////
//   // Check value routine  //
//   ///////////////////////
//   double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
//   double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
//   double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;

//   int thr_per_node = nthr/nnodes;

// #pragma omp parallel for
//   for (int thr=0; thr<nthr; thr++) {
//     int node = thr/thr_per_node;
//     multi_start[thr] = omp_get_wtime();
//     for (int i=0; i<num_vals; i++) {
//       double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end; 
//       double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end; 
//       double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end; 
//       eval_multi_UBspline_3d_d (multi_spline[node], x, y, z, multi_vals[thr]);
//     }
//     multi_end[thr] = omp_get_wtime();
//   }

// #pragma omp parallel for
//   for (int thr=0; thr<nthr; thr++) {
//     norm_start[thr] = omp_get_wtime();
//     for (int i=0; i<num_vals; i++) {
//       double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
//       double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
//       double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
//       for (int j=0; j<num_splines; j++)
// 	eval_UBspline_3d_d (norm_splines[j], x, y, z, &(norm_vals[thr][j]));
//     }
//     norm_end[thr] = omp_get_wtime();
//   }
  
//   double norm_avg=0.0, multi_avg=0.0;

//   for (int thr=0; thr<nthr; thr++) {
//     double norm_time   = (double)(norm_end[thr] - norm_start[thr] + rand_start - rand_end);
//     double multi_time  = (double)(multi_end[thr] - multi_start[thr] + rand_start - rand_end);
//     norm_avg += norm_time;
//     multi_avg += multi_time;
//   }
//   norm_avg  /= nthr;
//   multi_avg /= nthr;
//   double norm_speed  = (double) num_vals*num_splines / norm_avg;
//   double multi_speed = (double) num_vals*num_splines / multi_avg;

//   fprintf (stderr, "Normal value speed = %13.3f evaluations per second.\n", 
// 	   norm_speed);
//   fprintf (stderr, "Multi  value speed = %13.3f evaluations per second.\n", 
// 	   multi_speed);

  
//   ///////////////////////
//   // Check VGH routine //
//   ///////////////////////
//   #pragma omp parallel for
//   for (int thr=0; thr<nthr; thr++) {
//     int node = thr/thr_per_node;
//     multi_start[thr] = omp_get_wtime();
//     for (int i=0; i<num_vals; i++) {
//       double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
//       double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
//       double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
//       eval_multi_UBspline_3d_d_vgh 
// 	(multi_spline[node], x, y, z,  multi_vals[thr], 
// 	 multi_grads[thr], multi_hess[thr]);
//     }
//     multi_end[thr] = omp_get_wtime();
//   }

// #pragma omp parallel for
//   for (int thr=0; thr<nthr; thr++) {
//     norm_start[thr] = omp_get_wtime();
//     for (int i=0; i<num_vals; i++) {
//       double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
//       double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
//       double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
//       for (int j=0; j<num_splines; j++)
// 	eval_UBspline_3d_d_vgh (norm_splines[j], x, y, z, &(norm_vals[thr][j]),
// 				&(norm_grads[thr][3*j]), &(norm_hess[thr][9*j]));
//     }
//     norm_end[thr] = omp_get_wtime();
//   }

//   norm_avg = multi_avg = 0.0;
//   for (int thr=0; thr<nthr; thr++) {
//     double norm_time   = (double)(norm_end[thr] - norm_start[thr] + rand_start - rand_end);
//     double multi_time  = (double)(multi_end[thr] - multi_start[thr] + rand_start - rand_end);
//     norm_avg += norm_time;
//     multi_avg += multi_time;
//   }
//   norm_avg  /= nthr;
//   multi_avg /= nthr;
//   norm_speed  = (double) num_vals*num_splines / norm_avg;
//   multi_speed = (double) num_vals*num_splines / multi_avg;

//   fprintf (stderr, "Normal VGH   speed = %13.3f evaluations per second.\n", 
// 	   norm_speed);
//   fprintf (stderr, "Multi  VGH   speed = %13.3f evaluations per second.\n", 
// 	   multi_speed);


//   destroy_Bspline (multi_spline);
//   for (int i=0; i<num_splines; i++)
//     destroy_Bspline(norm_splines[i]); 
}


#endif

void
time_3d_real_double_all()
{
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
  multi_UBspline_3d_d *multi_spline;
  
  // First, create multispline
  multi_spline = create_multi_UBspline_3d_d 
    (x_grid, y_grid, z_grid, xBC, yBC, zBC, num_splines);

  double data[Nx*Ny*Nz];
  // Now, create normal splines and set multispline data
  for (int i=0; i<num_splines; i++) {
    for (int j=0; j<Nx*Ny*Nz; j++)
      data[j] = (drand48()-0.5);
    norm_splines[i] = create_UBspline_3d_d 
      (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
    set_multi_UBspline_3d_d (multi_spline, i, data);
  }
  
  // Now, test random values
  int num_vals = 10000;
  double multi_vals[num_splines], norm_vals[num_splines];
  double multi_grads[3*num_splines], norm_grads[3*num_splines];
  double multi_lapl[num_splines], norm_lapl[num_splines];
  double multi_hess[9*num_splines], norm_hess[9*num_splines];

  double rand_start, rand_end, norm_start, norm_end, multi_start, multi_end;

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
  multi_start = get_time();
  double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
  double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
  double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;

  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end; 
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end; 
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end; 
    eval_multi_UBspline_3d_d (multi_spline, x, y, z, multi_vals);
  }
  multi_end = get_time();

  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_d (norm_splines[j], x, y, z, &(norm_vals[j]));
  }
  norm_end = get_time();
  
  double norm_time   = (double)(norm_end - norm_start + rand_start - rand_end) ;
  double multi_time  = (double)(multi_end - multi_start + rand_start - rand_end) ;
  double norm_speed  = (double) num_vals*num_splines / norm_time;
  double multi_speed = (double) num_vals*num_splines / multi_time;
  fprintf (stderr, "Normal value speed = %13.3f evaluations per second.\n", 
	   norm_speed);
  fprintf (stderr, "Multi  value speed = %13.3f evaluations per second.\n", 
	   multi_speed);
  

  ///////////////////////
  // Check VGH routine //
  ///////////////////////
  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    eval_multi_UBspline_3d_d_vgh (multi_spline, x, y, z, 
				  multi_vals, multi_grads, multi_hess);
  }
  multi_end = get_time();

  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_d_vgh (norm_splines[j], x, y, z, &(norm_vals[j]),
			  &(norm_grads[3*j]), &(norm_hess[9*j]));
  }
  norm_end = get_time();

  norm_time   = (double)(norm_end - norm_start + rand_start - rand_end) ;
  multi_time  = (double)(multi_end - multi_start + rand_start - rand_end) ;
  norm_speed  = (double) num_vals*num_splines / norm_time;
  multi_speed = (double) num_vals*num_splines / multi_time;
  fprintf (stderr, "Normal VGH   speed = %13.3f evaluations per second.\n", 
	   norm_speed);
  fprintf (stderr, "Multi  VGH   speed = %13.3f evaluations per second.\n", 
	   multi_speed);

  destroy_Bspline (multi_spline);
  for (int i=0; i<num_splines; i++)
    destroy_Bspline(norm_splines[i]); 
}






void
time_3d_complex_double_all()
{
  int Nx=37; int Ny=37; int Nz = 37;
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
  multi_UBspline_3d_z *multi_spline;
  
  // First, create multispline
  multi_spline = create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC,
					     num_splines);

  complex_double data[Nx*Ny*Nz];
  // Now, create normal splines and set multispline data
  for (int i=0; i<num_splines; i++) {
    for (int j=0; j<Nx*Ny*Nz; j++)
      data[j] = (drand48()-0.5);// + (drand48()-0.5)*1.0i;
    norm_splines[i] = create_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
    set_multi_UBspline_3d_z (multi_spline, i, data);
  }
  
  // Now, test random values
  int num_vals = 100000;
  complex_double multi_vals[num_splines], norm_vals[num_splines];
  complex_double multi_grads[3*num_splines], norm_grads[3*num_splines];
  complex_double multi_lapl[num_splines], norm_lapl[num_splines];
  complex_double multi_hess[9*num_splines], norm_hess[9*num_splines];

  double rand_start, rand_end, norm_start, norm_end, multi_start, multi_end;

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
  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    eval_multi_UBspline_3d_z (multi_spline, x, y, z, multi_vals);
  }
  multi_end = get_time();

  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_z (norm_splines[j], x, y, z, &(norm_vals[j]));
  }
  norm_end = get_time();
  
  double norm_time   = (double)(norm_end - norm_start + rand_start - rand_end) ;
  double multi_time  = (double)(multi_end - multi_start + rand_start - rand_end) ;
  double norm_speed  = (double) num_vals*num_splines / norm_time;
  double multi_speed = (double) num_vals*num_splines / multi_time;
  fprintf (stderr, "Normal value speed = %13.3f evaluations per second.\n", norm_speed);
  fprintf (stderr, "Multi  value speed = %13.3f evaluations per second.\n", multi_speed);
  
  ///////////////////////
  // Check VGL routine //
  ///////////////////////
  // eval_multi_UBspline_3d_z_vgl (multi_spline, x, y, z, 
  // 				multi_vals, multi_grads, multi_lapl);
  // for (int j=0; j<num_splines; j++)
  //   eval_UBspline_3d_z_vgl (norm_splines[j], x, y, z, &(norm_vals[j]),
  // 			    &(norm_grads[3*j]), &(norm_lapl[j]));
  // for (int j=0; j<num_splines; j++) {
  //   // Check value
  //   if (zdiff(norm_vals[j], multi_vals[j], 1.0e-12))
  //     return -3;
    
  //   // Check gradients
  //   for (int n=0; n<3; n++) 
  //     if (zdiff (norm_grads[3*j+n], multi_grads[3*j+n], 1.0e-10))
  // 	return -4;
    
  //   // Check laplacian
  //   if (zdiff (norm_lapl[j], multi_lapl[j], 1.0e-10)) 
  //     return -5;
  // }


  ///////////////////////
  // Check VGH routine //
  ///////////////////////
  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    eval_multi_UBspline_3d_z_vgh (multi_spline, x, y, z, 
				  multi_vals, multi_grads, multi_hess);
  }
  multi_end = get_time();

  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_z_vgh (norm_splines[j], x, y, z, &(norm_vals[j]),
			  &(norm_grads[3*j]), &(norm_hess[9*j]));
  }
  norm_end = get_time();

  norm_time   = (double)(norm_end - norm_start + rand_start - rand_end) ;
  multi_time  = (double)(multi_end - multi_start + rand_start - rand_end) ;
  norm_speed  = (double) num_vals*num_splines / norm_time;
  multi_speed = (double) num_vals*num_splines / multi_time;
  fprintf (stderr, "Normal VGH   speed = %13.3f evaluations per second.\n", norm_speed);
  fprintf (stderr, "Multi  VGH   speed = %13.3f evaluations per second.\n", multi_speed);

  destroy_Bspline (multi_spline);
  for (int i=0; i<num_splines; i++)
    destroy_Bspline(norm_splines[i]);
 
}


void test_complex_double_vgh()
{
  int Nx=73; int Ny=91; int Nz = 29;
  int num_splines = 128;

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
  multi_UBspline_3d_z *multi_spline;
  
  // First, create multispline
  multi_spline = create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC,
					     num_splines);

  complex_double data[Nx*Ny*Nz];
  // Now, create normal splines and set multispline data
  for (int i=0; i<num_splines; i++) {
    for (int j=0; j<Nx*Ny*Nz; j++)
      data[j] = (drand48()-0.5);// + (drand48()-0.5)*1.0i;
    norm_splines[i] = create_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
    set_multi_UBspline_3d_z (multi_spline, i, data);
  }

  fprintf (stderr, "norm coef  = %1.14e + %1.14ei\n",
	   creal(norm_splines[19]->coefs[227]),
	   cimag(norm_splines[19]->coefs[227]));
  fprintf (stderr, "multi coef = %1.14e + %1.14ei\n",
	   creal(multi_spline->coefs[19+227*multi_spline->z_stride]),
	   cimag(multi_spline->coefs[19+227*multi_spline->z_stride]));
  
  // Now, test random values
  int num_vals = 100;
  complex_double multi_vals[num_splines], norm_vals[num_splines];
  complex_double multi_grads[3*num_splines], norm_grads[3*num_splines];
  complex_double multi_lapl[num_splines], norm_lapl[num_splines];
  complex_double multi_hess[9*num_splines], norm_hess[9*num_splines];
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    
    ///////////////////////
    // Check VGH routine //
    ///////////////////////
    eval_multi_UBspline_3d_z_vgh (multi_spline, x, y, z, 
				  multi_vals, multi_grads, multi_hess);
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_z_vgh (norm_splines[j], x, y, z, &(norm_vals[j]),
			  &(norm_grads[3*j]), &(norm_hess[9*j]));
    for (int j=0; j<num_splines; j++) {
      // Check value
      if (zdiff(norm_vals[j], multi_vals[j], 1.0e-12)) {
	fprintf (stderr, "Error!  norm_vals[j] = %1.14e + %1.14ei\n",
		 creal(norm_vals[j]), cimag(norm_vals[j]));
	fprintf (stderr, "       multi_vals[j] = %1.14e + %1.14ei\n",
		 creal(multi_vals[j]), cimag(multi_vals[j]));
      }
      // Check gradients
      for (int n=0; n<3; n++) {
	if (zdiff (norm_grads[3*j+n], multi_grads[3*j+n], 1.0e-12)) {
	  fprintf (stderr, "n=%d\n", n);
	  fprintf (stderr, "Error!  norm_grads[j] = %1.14e + %1.14ei\n",
		   creal(norm_grads[3*j+n]), cimag(norm_grads[3*j+n]));
	  fprintf (stderr, "       multi_grads[j] = %1.14e + %1.14ei\n",
		   creal(multi_grads[3*j+n]), cimag(multi_grads[3*j+n]));
	}
      }
      // Check hessian
      for (int n=0; n<9; n++) {
	if (zdiff (norm_hess[9*j+n], multi_hess[9*j+n], 1.0e-10)) {
	  fprintf (stderr, "Error!  norm_hess[j] = %1.14e + %1.14ei\n",
		   creal(norm_hess[9*j+n]), cimag(norm_hess[9*j+n]));
	  fprintf (stderr, "       multi_hess[j] = %1.14e + %1.14ei\n",
		   creal(multi_hess[9*j+n]), cimag(multi_hess[9*j+n]));
	}
      }
    }
  }

  num_vals = 100000;

  // Now do timing
  double norm_start, norm_end, multi_start, multi_end, rand_start, rand_end;
  rand_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
  }
  rand_end = get_time();

  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_z_vgh (norm_splines[j], x, y, z, &(norm_vals[j]),
			      &(norm_grads[3*j]), &(norm_hess[9*j]));
  }
  norm_end = get_time();

  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    eval_multi_UBspline_3d_z_vgh (multi_spline, x, y, z, multi_vals, multi_grads, multi_hess);
  }
  multi_end = get_time();

  fprintf (stderr, "Normal spline time = %1.5f\n",
	   (double)(norm_end-norm_start+rand_start-rand_end));
  fprintf (stderr, "Multi  spline time = %1.5f\n",
	   (double)(multi_end-multi_start+rand_start-rand_end));


}


void test_double()
{
  int Nx=73; int Ny=91; int Nz = 29;
  int num_splines = 201;

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
  multi_UBspline_3d_d *multi_spline;
  
  // First, create multispline
  multi_spline = create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC,
					     num_splines);
  
  double data[Nx*Ny*Nz];
  // Now, create normal splines and set multispline data
  for (int i=0; i<num_splines; i++) {
    for (int j=0; j<Nx*Ny*Nz; j++)
      data[j] = (drand48()-0.5);// + (drand48()-0.5)*1.0i;
    norm_splines[i] = create_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
    set_multi_UBspline_3d_d (multi_spline, i, data);
  }
  
  fprintf (stderr, "norm coef  = %1.14e\n",
	   norm_splines[19]->coefs[227]);
  fprintf (stderr, "multi coef = %1.14e\n",
	   multi_spline->coefs[19+227*multi_spline->z_stride]);
  
  // Now, test random values
  int num_vals = 100;
  double multi_vals[num_splines], norm_vals[num_splines];
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    
    eval_multi_UBspline_3d_d (multi_spline, x, y, z, 
			      multi_vals);
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_d (norm_splines[j], x, y, z, &(norm_vals[j]));
    for (int j=0; j<num_splines; j++) {
      // Check value
      double diff = norm_vals[j] - multi_vals[j];
      if (fabs(diff) > 1.0e-12) {
	fprintf (stderr, "Error!  norm_vals[j] = %1.14e\n",
		 norm_vals[j]);
	fprintf (stderr, "       multi_vals[j] = %1.14e\n",
		 multi_vals[j]);
      }
    }
  }
  
  num_vals = 100000;
  
  // Now do timing
  double norm_start, norm_end, multi_start, multi_end, rand_start, rand_end;
  rand_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
  }
  rand_end = get_time();
  
  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_d (norm_splines[j], x, y, z, &(norm_vals[j]));
  }
  norm_end = get_time();
  
  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    eval_multi_UBspline_3d_d (multi_spline, x, y, z, multi_vals);
  }
  multi_end = get_time();
  
  fprintf (stderr, "Normal spline time = %1.5f\n",
	   (double)(norm_end-norm_start+rand_start-rand_end));
  fprintf (stderr, "Multi  spline time = %1.5f\n",
	   (double)(multi_end-multi_start+rand_start-rand_end));
  
}



void test_double_vgh()
{
  int Nx=73; int Ny=91; int Nz = 29;
  int num_splines = 128;

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
  multi_UBspline_3d_d *multi_spline;
  
  // First, create multispline
  multi_spline = create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC,
					     num_splines);
  
  double data[Nx*Ny*Nz];
  // Now, create normal splines and set multispline data
  for (int i=0; i<num_splines; i++) {
    for (int j=0; j<Nx*Ny*Nz; j++)
      data[j] = (drand48()-0.5);// + (drand48()-0.5)*1.0i;
    norm_splines[i] = create_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
    set_multi_UBspline_3d_d (multi_spline, i, data);
  }
  
  fprintf (stderr, "norm coef  = %1.14e\n",
	   norm_splines[19]->coefs[227]);
  fprintf (stderr, "multi coef = %1.14e\n",
	   multi_spline->coefs[19+227*multi_spline->z_stride]);
  
  // Now, test random values
  int num_vals = 100;
  double multi_vals[num_splines], norm_vals[num_splines];
  double multi_grads[3*num_splines], norm_grads[3*num_splines];
  double multi_hess[9*num_splines], norm_hess[9*num_splines];
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    
    eval_multi_UBspline_3d_d_vgh (multi_spline, x, y, z, 
				  multi_vals, multi_grads, multi_hess);
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_d_vgh (norm_splines[j], x, y, z, &(norm_vals[j]),
			      &(norm_grads[3*j]), &(norm_hess[9*j]));
    for (int j=0; j<num_splines; j++) {
      // Check value
      double diff = norm_vals[j] - multi_vals[j];
      if (fabs(diff) > 1.0e-12) {
	fprintf (stderr, "j = %d\n", j);
	fprintf (stderr, "Error!  norm_vals[j] = %1.14e\n",
		 norm_vals[j]);
	fprintf (stderr, "       multi_vals[j] = %1.14e\n",
		 multi_vals[j]);
      }
      // Check gradients
      for (int n=0; n<3; n++) {
	diff = norm_grads[3*j+n] - multi_grads[3*j+n];
	if (fabs(diff) > 1.0e-12) {
	  fprintf (stderr, "n=%d\n", n);
	  fprintf (stderr, "Error!  norm_grads[j] = %1.14e\n",
		   norm_grads[3*j+n]);
	  fprintf (stderr, "       multi_grads[j] = %1.14e\n",
		   multi_grads[3*j+n]);
	}
      }
      // Check hessian
      for (int n=0; n<9; n++) {
	diff = norm_hess[9*j+n] - multi_hess[9*j+n];
	if (fabs(diff) > 1.0e-10) {
	  fprintf (stderr, "Error!  norm_hess[j] = %1.14e\n",
		   norm_hess[9*j+n]);
	  fprintf (stderr, "       multi_hess[j] = %1.14e\n",
		   multi_hess[9*j+n]);
	}
      }
    }
  }
  
  num_vals = 100000;
  
  // Now do timing
  double norm_start, norm_end, multi_start, multi_end, rand_start, rand_end;
  rand_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
  }
  rand_end = get_time();
  
  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    
    for (int j=0; j<num_splines; j++)
      eval_UBspline_3d_d_vgh (norm_splines[j], x, y, z, &(norm_vals[j]),
			      &(norm_grads[3*j]), &(norm_hess[9*j]));
  }
  norm_end = get_time();
  
  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  double x = rx*x_grid.start + (1.0-rx)*x_grid.end;
    double ry = drand48();  double y = ry*y_grid.start + (1.0-ry)*y_grid.end;
    double rz = drand48();  double z = rz*z_grid.start + (1.0-rz)*z_grid.end;
    eval_multi_UBspline_3d_d_vgh (multi_spline, x, y, z, multi_vals, multi_grads, multi_hess);
  }
  multi_end = get_time();
  
  fprintf (stderr, "Normal spline time = %1.5f\n",
	   (double)(norm_end-norm_start+rand_start-rand_end));
  fprintf (stderr, "Multi  spline time = %1.5f\n",
	   (double)(multi_end-multi_start+rand_start-rand_end));
  
}



int 
time_1d_NUB_complex_double_all()
{
  int Nx=100;
  int num_splines = 128*36;

  NUgrid *x_grid = create_log_grid (1.0e-4, 3.0, Nx);
  //  for (int i=0; i<Nx; i++) 
  //  fprintf (stderr, "%1.8e\n", x_grid->points[i]);

  BCtype_z xBC;
  // xBC.lCode = xBC.rCode = NATURAL;
  xBC.lCode = DERIV1; xBC.lVal_r = 2.3; xBC.lVal_i = 1.1;
  xBC.rCode = DERIV1; xBC.rVal_r = -2.3; xBC.rVal_i = -1.1;
  

  // First, create splines the normal way
  NUBspline_1d_z* norm_splines[num_splines];
  multi_NUBspline_1d_z *multi_spline;
  
  // First, create multispline
  multi_spline = create_multi_NUBspline_1d_z (x_grid, xBC, num_splines);

  complex_double data[Nx];
  // Now, create normal splines and set multispline data
  for (int i=0; i<num_splines; i++) {
    for (int j=0; j<Nx; j++)
      data[j] = (drand48()-0.5);// + (drand48()-0.5)*1.0i;

    xBC.lVal_r = drand48(); xBC.lVal_i = drand48();
    xBC.rVal_r = drand48(); xBC.rVal_i = drand48();

    norm_splines[i] = create_NUBspline_1d_z (x_grid, xBC, data);
    //set_multi_NUBspline_1d_z (multi_spline, i, data);
    set_multi_NUBspline_1d_z_BC (multi_spline, i, data, xBC);
  }
  
  // Now, test random values
  int num_vals = 100000;
  complex_double  multi_vals[num_splines], norm_vals [num_splines];
  complex_double multi_grads[num_splines], norm_grads[num_splines];
  complex_double  multi_lapl[num_splines], norm_lapl [num_splines];

  double multi_start, multi_end, norm_start, norm_end;
  

  //////////////////////////
  // Time value routine   //
  //////////////////////////
  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  
    double x = rx*x_grid->start + (1.0-rx)*x_grid->end;

    eval_multi_NUBspline_1d_z (multi_spline, x, multi_vals);
  }
  multi_end = get_time();

  norm_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  
    double x = rx*x_grid->start + (1.0-rx)*x_grid->end;

    for (int j=0; j<num_splines; j++)
      eval_NUBspline_1d_z (norm_splines[j], x, &(norm_vals[j]));
  }
  norm_end = get_time();
  double dt = (double)(multi_end - multi_start) ;
  double multi_speed = (double)num_vals * (double)num_splines/ dt; 
  fprintf (stderr, "1D complex nonuniform multi-spline speed = %9.2f\n",
	   multi_speed);


  //////////////////////////
  // Time VGL routine   //
  //////////////////////////
  multi_start = get_time();
  for (int i=0; i<num_vals; i++) {
    double rx = drand48();  
    double x = rx*x_grid->start + (1.0-rx)*x_grid->end;
    eval_multi_NUBspline_1d_z_vgl (multi_spline, x, multi_vals, multi_grads, multi_lapl);
  }
  multi_end = get_time();

  /* norm_start = get_time(); */
  /* for (int i=0; i<num_vals; i++) { */
  /*   double rx = drand48();   */
  /*   double x = rx*x_grid->start + (1.0-rx)*x_grid->end; */

  /*   for (int j=0; j<num_splines; j++) */
  /*     eval_NUBspline_1d_z (norm_splines[j], x, &(norm_vals[j])); */
  /* } */
  /* norm_end = get_time(); */
  dt = (double)(multi_end - multi_start) ;
  multi_speed = (double)num_vals * (double)num_splines/ dt; 
  fprintf (stderr, "1D complex nonuniform multi-spline speed = %9.2f\n",
	   multi_speed);


 return 0;
}




void PrintPassFail (int code)
{
  char green[100], normal[100], red[100];
  snprintf (green, 100,  "%c[0;32;47m", 0x1b);
  snprintf (normal, 100, "%c[0;30;47m", 0x1b);
  snprintf (red,    100, "%c[0;31;47m", 0x1b);

  if (code == 0) 
    fprintf (stderr, "PASSED\n");
  else 
    fprintf (stderr, "FAILED:  code = %d\n", code);
}


int main()
{
  // time_1d_NUB_complex_double_all();
// #ifdef _OPENMP
//   fprintf (stderr, "Timing 3D double-precision evaluation speed with OpenMP:\n");
//   time_3d_real_double_omp();
// #endif
//   fprintf (stderr, "Timing 3D complex single-precision evaluation speed:\n"); 
//   time_3d_complex_float_all(); 
//   fprintf (stderr, "Timing 3D single-precision evaluation speed:\n");
//   time_3d_real_float_all();
  fprintf (stderr, "Timing 3D double-precision evaluation speed:\n");
  time_3d_real_double_all();
  fprintf (stderr, "Timing 3D complex double-precision evaluation speed:\n");
  time_3d_complex_double_all();
//  test_3d_double_all();
}
