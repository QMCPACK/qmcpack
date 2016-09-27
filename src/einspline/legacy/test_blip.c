//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign 
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#include "blip_create.h"
#include "bspline.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

inline
double dot (double a[3], double b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline
void cross (double a[3], double b[3], double axb[3])
{
  axb[0] = a[1]*b[2]-a[2]*b[1];
  axb[1] = a[2]*b[0]-a[0]*b[2];
  axb[2] = a[0]*b[1]-a[1]*b[0];
}

inline 
void recip_lattice (double lattice[9], double recip[9])
{
  double *a0 = lattice+0; double *b0 = recip+0;
  double *a1 = lattice+3; double *b1 = recip+3;
  double *a2 = lattice+6; double *b2 = recip+6;
  double a1xa2[3];
  cross (a1, a2, a1xa2);
  double vol = dot (a0, a1xa2);
  double prefactor = 2.0 * M_PI / vol;
  
  cross (a1, a2, b0);
  cross (a2, a0, b1);
  cross (a0, a1, b2);
  for (int i=0; i<3; i++) {
    b0[i] *= prefactor;
    b1[i] *= prefactor;
    b2[i] *= prefactor;
  }
}


void
test_blip_1d_s (double a, double Gcut, double ratio)
{
  double lattice[9] = { 0.0,  a ,  a  ,
                         a , 0.0,  a  ,
                         a ,  a , 0.0 };
  double  recip[9];
  recip_lattice (lattice, recip);
  double *b0 = recip+0;
  double *b1 = recip+3;
  double *b2 = recip+6;
  int x_max = (int)ceil(Gcut / sqrt (dot (b0, b0)));
  int y_max = (int)ceil(Gcut / sqrt (dot (b1, b1)));
  int z_max = (int)ceil(Gcut / sqrt (dot (b2, b2)));

  int numG = 0;
  double G[3], G0[3], G1[3], G2[3];
  // Count G-vectors
  for (int ix=-x_max; ix<=x_max; ix++) {
    G0[0] = (double)ix * b0[0];
    G0[1] = (double)ix * b0[1];
    G0[2] = (double)ix * b0[2];
    for (int iy=-y_max; iy<=y_max; iy++) {
      G1[0] = (double)iy * b1[0];
      G1[1] = (double)iy * b1[1];
      G1[2] = (double)iy * b1[2];
      for (int iz=-z_max; iz<=z_max; iz++) {
	G2[0] = (double)iz * b2[0];
	G2[1] = (double)iz * b2[1];
	G2[2] = (double)iz * b2[2];
	
	G[0] = G0[0] + G1[0] + G2[0];
	G[1] = G0[1] + G1[1] + G2[1];
	G[2] = G0[2] + G1[2] + G2[2];
	
	double gmag = dot (G, G);
	if (gmag < Gcut*Gcut)
	  numG++;
      }
    }
  }
  fprintf (stderr, "There are %d G-vectors\n", numG);
  double *Gvecs = malloc (sizeof(double)*numG*3);
  complex_float *coefs = malloc (sizeof(complex_float)*numG);
  
  numG = 0;
  for (int ix=-x_max; ix<=x_max; ix++) {
    G0[0] = (double)ix * b0[0];
    G0[1] = (double)ix * b0[1];
    G0[2] = (double)ix * b0[2];
    for (int iy=-y_max; iy<=y_max; iy++) {
      G1[0] = (double)iy * b1[0];
      G1[1] = (double)iy * b1[1];
      G1[2] = (double)iy * b1[2];
      for (int iz=-z_max; iz<=z_max; iz++) {
	G2[0] = (double)iz * b2[0];
	G2[1] = (double)iz * b2[1];
	G2[2] = (double)iz * b2[2];
	
	G[0] = G0[0] + G1[0] + G2[0];
	G[1] = G0[1] + G1[1] + G2[1];
	G[2] = G0[2] + G1[2] + G2[2];
	
	double gmag = dot (G, G);
	if (gmag < Gcut*Gcut) {
	  Gvecs[numG*3+0] = G[0];
	  Gvecs[numG*3+1] = G[1];
	  Gvecs[numG*3+2] = G[2];
	  coefs[numG] = (float)drand48() + (float)drand48()*1.0fi;
	  numG++;
	}
      }
    }
  }

  UBspline_3d_s *blip, *interp;
  blip   = create_blip_3d_s (lattice, Gvecs, coefs, numG, 4.0, true);

  // Now, evaluate spline on a line through the box
  FILE *fout = fopen ("blip_1d_s.dat", "w");
  double y = 0.2; 
  double z = 0.8;
  double u[3], r[3];
  u[1] = y;
  u[2] = z;
  for (double x=0.0; x<1.0; x+=0.001) {
    u[0] = x;
    r[0] = u[0]*lattice[0] + u[1]*lattice[3] + u[2]*lattice[6];
    r[1] = u[0]*lattice[1] + u[1]*lattice[4] + u[2]*lattice[7];
    r[2] = u[0]*lattice[2] + u[1]*lattice[5] + u[2]*lattice[8];

    float val, derivs[3];
    eval_UBspline_3d_s_vg (blip, x, 0.2, 0.8, &val, derivs);
    
    double sum = 0.0;
    double gsum[3] = { 0.0, 0.0, 0.0}; 
    double gex[3];
    gsum[0] = 0.0; gsum[1] = 0.0; gsum[2] = 0.0;
    for (int i=0; i<numG; i++) {
      double phase = dot (Gvecs+3*i, r);
      complex_float e2iGr = (float)cos(phase) + 1.0fi*sin(phase);
      sum += crealf (coefs[i] * e2iGr);
      gsum[0] += crealf (coefs[i] *e2iGr *1.0fi * Gvecs[3*i+0]);
      gsum[1] += crealf (coefs[i] *e2iGr *1.0fi * Gvecs[3*i+1]);
      gsum[2] += crealf (coefs[i] *e2iGr *1.0fi * Gvecs[3*i+2]);
    }
    float grad[3];
    grad[0] = recip[0]*derivs[0] + recip[1]*derivs[1] + recip[2]*derivs[2];
    grad[1] = recip[3]*derivs[0] + recip[4]*derivs[1] + recip[5]*derivs[2];
    grad[2] = recip[6]*derivs[0] + recip[7]*derivs[1] + recip[8]*derivs[2];
    grad[0]/=(2.0*M_PI); grad[1]/=(2.0*M_PI); grad[2]/=(2.0*M_PI);

    fprintf (fout, "%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n", x, val, sum, grad[0], grad[1], grad[2],
	     gsum[0], gsum[1], gsum[2]);
  }
  fclose (fout);

  destroy_Bspline (blip);
  free (Gvecs);
  free (coefs);
}


int main()
{
  test_blip_1d_s(3.0, 20.0, 1.0);
}
