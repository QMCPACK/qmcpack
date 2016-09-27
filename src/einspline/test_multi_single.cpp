//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#include "multi_bspline.h"
#include "multi_nubspline.h"
#include "bspline.h"
#include "nubspline.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h>
#include <getopt.h>

inline double get_time()
{
  return omp_get_wtime();
}


/** A simplified SPOSet using einspline in single precision
 */
struct EinsplineSet
{
  int Nx, Ny, Nz;
  int num_splines;
  ///spline engine
  multi_UBspline_3d_s *multi_spline;

  EinsplineSet(int nx, int ny, int nz, int ns, bool init)
    : Nx(nx),Ny(ny),Nz(nz),num_splines(ns)
  {
    Ugrid x_grid, y_grid, z_grid;
    x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num = Nx;
    y_grid.start = 0.0; y_grid.end = 1.0; y_grid.num = Ny;
    z_grid.start = 0.0; z_grid.end = 1.0; z_grid.num = Nz;

    BCtype_s xBC, yBC, zBC;
    xBC.lCode = xBC.rCode = PERIODIC;
    yBC.lCode = yBC.rCode = PERIODIC;
    zBC.lCode = zBC.rCode = PERIODIC;

    // First, create multispline
    multi_spline = create_multi_UBspline_3d_s(x_grid, y_grid, z_grid, xBC, yBC, zBC, num_splines);

    if(init)
    {
      float data[Nx*Ny*Nz];
      // Now, create normal splines and set multispline data
      for (int i=0; i<num_splines; i++) 
      {
        for (int j=0; j<Nx*Ny*Nz; j++)
          data[j] = (drand48()-0.5);// + (drand48()-0.5)*1.0i;
        set_multi_UBspline_3d_s (multi_spline, i, data);
      }
    }
    else
    {
      std::fill(multi_spline->coefs,multi_spline->coefs+multi_spline->coefs_size,0.0);
    }
  }

  ~EinsplineSet()
  {
    free(multi_spline);
  }


  inline void evaluate_v(float x, float y, float z, float* multi_vals) const
  {
    eval_multi_UBspline_3d_s (multi_spline, x, y, z, multi_vals);
  }

  inline void evaluate_vgh(float x, float y, float z, float* restrict multi_vals, float* restrict multi_g, float* restrict multi_h) const
  {
    eval_multi_UBspline_3d_s_vgh(multi_spline, x, y, z, multi_vals, multi_g, multi_h);
  }

};


template<typename T>
inline void randomize(T* pos, int n)
{
  for(int i=0; i<n*3; ++i) pos[i]=drand48();
}

int main(int argc, char **argv)
{
  int nx=32; 
  int ny=32; 
  int nz=32;
  int num_splines = 128;
  int niters=10;
  int nparticles =num_splines;

  int opt;
  while((opt = getopt(argc, argv, "hg:x:y:z:i:s:p:")) != -1)
  {
    switch(opt)
    {
    case 'h':
      printf("[-g grid| -x grid_x -y grid_y -z grid_z] -s states -p particles -i iterations  \n");
      return 1;
    case 'g': //set the grid a cubic box
      nx=ny=nz=atoi(optarg);
      break;
    case 'x'://set the xgrid
      nx=atoi(optarg);
      break;
    case 'y'://set the ygrid
      ny=atoi(optarg);
      break;
    case 'z'://set the zgrid
      nz=atoi(optarg);
      break;
    case 's': //number of splines
      num_splines=atoi(optarg);
      break;
    case 'p': //number of particles
      nparticles=atoi(optarg);
      break;
    case 'i': //number of iterations
      niters=atoi(optarg);
      break;
    }
  }

  //if true, initialize random values, debugging purpose
  bool init_random=false;
  EinsplineSet orb(nx,ny,nz,num_splines,init_random);

#pragma omp parallel
  {
    float pos[3*nparticles], sphere[36];
    float vals[num_splines];
    float grads[3*num_splines];
    float hess[9*num_splines];
    float fx=2.0/static_cast<float>(nx);
    float fy=2.0/static_cast<float>(ny);
    float fz=2.0/static_cast<float>(nz);


    for(int iter=0; iter<niters; ++iter)
    {
      randomize(pos,nparticles);

      //stage 1: diffusion
      //timer1.start();
      for(int iat=0,i=0;iat<nparticles; ++iat,i+=3) 
        orb.evaluate_vgh(pos[i],pos[i+1],pos[i+2],vals,grads,hess);
      //timer1.end();

      //sleep(drand48());

      //stage 2: PP evaluations 
      //timer1.start();
      for(int iat=0,i=0;iat<nparticles; ++iat,i+=3) 
      {
        randomize(sphere,12);
        for(int k=0,kk=0; k<12; ++k,kk+=3)
          orb.evaluate_v(pos[i]+fx*sphere[kk],pos[i+1]+fy*sphere[kk+1],pos[i+2]+fz*sphere[kk+2],vals);
      }
      //timer1.end();
    }
  }

  return 0;
}
