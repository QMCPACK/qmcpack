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

#include <string.h>
#include <stdio.h>
#include "multi_bspline_create.h"
#include "multi_bspline_copy.h"
#include "multi_bspline_eval_d.h"
#include "bspline_create.h"

#ifdef __cplusplus
extern "C" {
#endif

  void copy_UBspline_3d_d(multi_UBspline_3d_d* multi, int i
      , const UBspline_3d_d* single, const int* offset, const int* N)
  {
    intptr_t z_stride=multi->z_stride;

    for(int ix=0; ix<N[0]; ++ix)
      for(int iy=0; iy<N[1]; ++iy)
      {
        intptr_t out=ix*multi->x_stride+iy*multi->y_stride+i;
        intptr_t in =(ix+offset[0])*single->x_stride+(iy+offset[1])*single->y_stride+offset[2];
        for(int iz=0; iz<N[2]; ++iz, in++, out +=z_stride)
        {
          multi->coefs[out]=single->coefs[in];
        }
      }
  }

  void copy_UBspline_3d_d_s(multi_UBspline_3d_s* multi, int i
      , const UBspline_3d_d* single, const int* offset, const int* N)
  {
    intptr_t x_stride_in=single->x_stride;
    intptr_t y_stride_in=single->y_stride;
    intptr_t x_stride_out=multi->x_stride;
    intptr_t y_stride_out=multi->y_stride;
    intptr_t z_stride_out=multi->z_stride;
    intptr_t offset0=(intptr_t)offset[0];
    intptr_t offset1=(intptr_t)offset[1];
    intptr_t offset2=(intptr_t)offset[2];
    const intptr_t istart=(intptr_t)i;
    const intptr_t n0=N[0],n1=N[1],n2=N[2];
    for(intptr_t ix=0; ix<n0; ++ix)
      for(intptr_t iy=0; iy<n1; ++iy)
      {
        float* restrict out=multi->coefs+ix*x_stride_out+iy*y_stride_out+istart;
        const double* restrict in =single->coefs+(ix+offset0)*x_stride_in+(iy+offset1)*y_stride_in+offset2;
        for(intptr_t iz=0; iz<n2; ++iz)
        {
          out[iz*z_stride_out]=(float)in[iz];
        }
      }
//
//    for(int ix=0; ix<N[0]; ++ix)
//      for(int iy=0; iy<N[1]; ++iy)
//      {
//        intptr_t out=ix*multi->x_stride+iy*multi->y_stride+i;
//        intptr_t in =(ix+offset[0])*single->x_stride+(iy+offset[1])*single->y_stride+offset[2];
//        for(int iz=0; iz<N[2]; ++iz, in++, out +=z_stride)
//        {
//          multi->coefs[out]=(float)single->coefs[in];
//        }
//      }
  }

  // Create 3D uniform single-precision, real Bspline
  multi_UBspline_3d_s *
    copy_multi_UBspline_3d_s (multi_UBspline_3d_s* spline)
    {
      multi_UBspline_3d_s *clone=create_multi_UBspline_3d_s(
          spline->x_grid, spline->y_grid, spline->z_grid
          ,spline->xBC, spline->yBC, spline->zBC
          ,spline->num_splines);
      memcpy(clone->coefs,spline->coefs,sizeof(float)*clone->coefs_size);
      return clone;
    }

  // just checking not used
  multi_UBspline_3d_d *
    copy_multi_UBspline_3d_d (multi_UBspline_3d_d* spline)
    {
      fprintf (stderr, "Who is using copy_multi_UBspline_3d_s?\n");
      abort();
      return NULL;
//      //first copy the grid
//      Ugrid x_grid=spline->x_grid;
//      Ugrid y_grid=spline->y_grid;
//      Ugrid z_grid=spline->z_grid;
//      
//      int ngy_i=(int)(0.21*y_grid.delta_inv);
//      int ngy_f=(int)(0.46*y_grid.delta_inv)+1;
//      y_grid.start=(double)(ngy_i)*y_grid.delta;
//      y_grid.end=(double)(ngy_f)*y_grid.delta;
//      y_grid.num=ngy_f-ngy_i;
//
//      int ngz_i=(int)(0.2*z_grid.delta_inv);
//      int ngz_f=(int)(0.5*z_grid.delta_inv)+1;
//      z_grid.start=(double)(ngz_i)*z_grid.delta;
//      z_grid.end=(double)(ngz_f)*z_grid.delta;
//      z_grid.num=ngz_f-ngz_i;
//
//      //check the boundary condition
//      int Nx_f=x_grid.num+3;
//      int Ny_f=y_grid.num+3;
//      int Nz_f=z_grid.num+3;
//
//      printf("start=%16.8f end=%16.8f delta=%16.8f ngz=%d\n"
//          ,spline->z_grid.start, spline->z_grid.end, spline->z_grid.delta, spline->z_grid.num);
//
//      printf("start=%16.8f end=%16.8f delta=%16.8f ngz=%d\n"
//          ,z_grid.start, z_grid.end, z_grid.delta, z_grid.num);
//
//      multi_UBspline_3d_d *clone=create_multi_UBspline_3d_d(
//          spline->x_grid, spline->y_grid, spline->z_grid
//          ,spline->xBC, spline->yBC, spline->zBC
//          ,spline->num_splines);
//
//      multi_UBspline_3d_d *small=create_multi_UBspline_3d_d(
//          x_grid, y_grid, z_grid
//          ,spline->xBC, spline->yBC, spline->zBC
//          ,spline->num_splines);
//
//
//      for(int ix=0; ix<Nx_f; ++ix)
//        for(int iy=0; iy<Ny_f; ++iy)
//          for(int iz=0; iz<Nz_f; ++iz)
//          {
//            intptr_t p_out=ix*small->x_stride+iy*small->y_stride+iz*small->z_stride;
//            intptr_t p_in =ix*spline->x_stride+(iy+ngy_i)*spline->y_stride+(iz+ngz_i)*spline->z_stride;
//            for(int i=0; i<spline->num_splines; ++i) small->coefs[p_out++]=spline->coefs[p_in++];
////
////            small->coeffs[ix*small->x_stride+iy*small->y_stride+iz]=
////            spline->coeffs[ix*spline->x_stride+iy*spline->y_stride+iz];
////            intptr_t offset = small->iy*Mz+iz;
//          }
//
//
//      memcpy(clone->coefs,spline->coefs,sizeof(double)*clone->coefs_size);
//
//      int num=spline->num_splines;
//      double vorg[num];
//      double vnew[num];
//
//      double x=0.001;
//      int ok=1;
//      int ntest=0;
//      while(x<1.0)
//      {
//        double y=y_grid.start+0.001;
//        while(y<y_grid.end)
//        {
//          double z=z_grid.start+1e-12;
//          while(z<z_grid.end)
//          {
//            eval_multi_UBspline_3d_d(spline,x,y,0.3,vorg);
//            eval_multi_UBspline_3d_d(small,x,y,0.3,vnew);
//
//            for(int i=0; i<num; ++i)
//              if(fabs(vorg[i]-vnew[i])>1e-12)
//              {
//                printf("%d %16.12f %16.12f %16.12f\n",i,vorg[i],vnew[i],vorg[i]-vnew[i]);
//                ok=0;
//              }
//            //  printf("%d %16.12f %16.12f %16.12f\n",i,vorg[i],vnew[i],vorg[i]-vnew[i]);
//            z+=0.02;
//            ntest++;
//          }
//          y+=0.07;
//        }
//        x += 0.03;
//      }
//
//      if(ok) 
//        printf("Everything is good with %d points\n",ntest);
//      else
//        printf("Something is wrong\n");
//      free(small);
//      return clone;
    }

  void copy_UBspline_1d_d(multi_UBspline_1d_d* multi, int i
      , const UBspline_1d_d* single, const int offset, const int N)
  {
    //fprintf(stdout,"debug xstride %ld i %d N %d \n", multi->x_stride, i, offset);
    for(int ix=0; ix<N; ++ix)
    {
      intptr_t out=ix*multi->x_stride+i;
      intptr_t in =ix+offset;
      multi->coefs[out]=single->coefs[in];
    }
  }

  void copy_UBspline_1d_d_s(multi_UBspline_1d_s* multi, int i
      , const UBspline_1d_d* single, const int offset, const int N)
  {
    for(int ix=0; ix<N; ++ix)
    {
      intptr_t out=ix*multi->x_stride+i;
      intptr_t in =ix+offset;
      multi->coefs[out]=(float)single->coefs[in];
    }
  }

#ifdef __cplusplus
}
#endif
