//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//		      Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign 
//
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////



/*
 * Copyright (C) 2004 PWSCF group 
 * Copyright (C) 2007 QMCPACK developers
 *
 * @author Jeongnim Kim http://www.mcc.uiuc.edu/qmcpack/
 * @brief Implements generic hdf5 interfaces for plane wave codes and qmcpack
 *
 * - pwhdf_open_file: open hdf5 file
 * - pwhdf_close_file : close hdf5 file
 * - pwhdf_open_eigg : open eigenstates
 * - pwhdf_close_eigg : close eigenstates
 * - pwhdf_open_eigr : open eigenstates_nx_ny_nz
 * - pwhdf_close_eigr : close eigenstates_nx_ny_nz
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "c_defs.h"
#include "hdf5.h"

/* file handler */
static hid_t h_file=-1;
/* main handler */
static hid_t h_main=-1;
/* twist angle handler */
static hid_t h_twist=-1;
/* number of fft grid */
static int h_ngrid[3];
/* number of real-space grids */
static int h_ngridtot=0;
/* check for gamma */
static int is_gamma=0;
void F77_FUNC_(pwhdf_open_file,PWHDF_OPEN_FILE)(const char* fname, const int* length)
{
  char * hfname = ( char * ) malloc( (*length) + 1 ) ;
  memcpy( hfname , fname , *length ) ;
  hfname[*length] = '\0' ; 

  if(h_file>=0) H5Fclose(h_file);

  h_file = H5Fcreate(hfname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  /* impelements version 1.00 hdf5 format */
  int version[]={1,10};
  hsize_t dim=2;
  hid_t dataspace= H5Screate_simple(1, &dim, NULL);
  hid_t dataset= H5Dcreate(h_file, "version", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,version);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  free(hfname);
}
void F77_FUNC_(pwhdf_close_file,PWHDF_CLOSE_FILE)()
{
  if(h_file>=0) H5Fclose(h_file);
  h_file=-1;
}

/** write basisset: number of plane waves, plane wave coefficients
 */
void F77_FUNC_(pwhdf_write_basis,PWHDF_WRITE_BASIS)(const int* ig, 
    const double* gcart, const int* ngtot)
{
  int ng=*ngtot;

  hid_t h_basis = H5Gcreate(h_file,"basis",0);
  hsize_t dim=1;
  hid_t dataspace= H5Screate_simple(1, &dim, NULL);
  hid_t dataset= H5Dcreate(h_basis, "num_planewaves", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ngtot);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  hsize_t dims[2];
  dims[0] = ng;
  dims[1] = 3;
  dataspace  = H5Screate_simple(2, dims, NULL);

  dataset =  H5Dcreate(h_basis, "multipliers", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ig);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  dataspace  = H5Screate_simple(2, dims, NULL);
  dataset =  H5Dcreate(h_basis, "planewaves", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,gcart);

  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Gclose(h_basis);
}


/** write basisset: number of plane waves, plane wave coefficients
void F77_FUNC_(pwhdf_write_basis,PWHDF_WRITE_BASIS)(const double* g, const int* igtog, const int* ngtot)
{
  int ng=*ngtot;
  int *ig=(int*)malloc(3*ng*sizeof(int));
  for(int i=0,i3=0; i<ng; i++)
  {
    int cur=3*(igtog[i]-1);
    ig[i3++]=(int)g[cur++];
    ig[i3++]=(int)g[cur++];
    ig[i3++]=(int)g[cur++];
  }

  hid_t h_basis = H5Gcreate(h_file,"basis",0);
  hsize_t dim=1;
  hid_t dataspace= H5Screate_simple(1, &dim, NULL);
  hid_t dataset= H5Dcreate(h_basis, "num_planewaves", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ngtot);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  hsize_t dims[2];
  dims[0] = ng;
  dims[1] = 3;
  dataspace  = H5Screate_simple(2, dims, NULL);
  dataset =  H5Dcreate(h_basis, "planewaves", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ig);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  H5Gclose(h_basis);

  free(ig);
}
 */
void F77_FUNC_(pwhdf_write_parameters,PWHDF_WRITE_PARAMETERS)(
    const int* nelec, const int* nspin, const int* nband, const int* nk,
    const double* ecut, const double* alat, const double* at)
{

  hid_t h_param = H5Gcreate(h_file,"parameters",0);
  hsize_t dim=1;
  hid_t dataspace= H5Screate_simple(1, &dim, NULL);
  hid_t dataset= H5Dcreate(h_param, "num_spins", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,nspin);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  dataspace= H5Screate_simple(1, &dim, NULL);
  dataset= H5Dcreate(h_param, "num_electrons", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,nelec);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  dataspace= H5Screate_simple(1, &dim, NULL);
  dataset= H5Dcreate(h_param, "num_bands", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,nband);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  dataspace= H5Screate_simple(1, &dim, NULL);
  dataset= H5Dcreate(h_param, "num_twists", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,nk);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  int iscomplex=1;
  dataspace= H5Screate_simple(1, &dim, NULL);
  dataset= H5Dcreate(h_param, "complex_coefficients", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&iscomplex);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  dataspace= H5Screate_simple(1, &dim, NULL);
  dataset= H5Dcreate(h_param, "maximum_ecut", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ecut);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  double lattice[9];
  for(int i=0; i<9; i++) lattice[i]=(*alat)*at[i];
  hsize_t dims[]={3,3};
  dataspace  = H5Screate_simple(2, dims, NULL);
  dataset =  H5Dcreate(h_param, "lattice", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,lattice);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  H5Gclose(h_param);
}

/* open mainbody:eigenstates */
void F77_FUNC_(pwhdf_open_eigg,PWHDF_OPEN_EIGG)()
{
  if(h_main>=0) H5Gclose(h_main);
  h_main = H5Gcreate(h_file,"eigenstates",0);
}

/* close eigenstates */
void F77_FUNC_(pwhdf_close_eigg,PWHDF_CLOSE_EIGG)()
{
  if(h_main>=0) H5Gclose(h_main);
  h_main=-1;
}

/* open twist# */
void F77_FUNC_(pwhdf_open_twist,PWHDF_OPEN_TWIST)(const int* ik, const double *xk,
    const int* nband, const int* nspin)
{
  char twistname[16];
  sprintf(twistname,"twist%i",(*ik)-1);
  if(h_twist>=0) H5Gclose(h_twist);
  h_twist = H5Gcreate(h_main,twistname,0);

  /* write twist_angle */
  hsize_t dim=3;
  hid_t dataspace= H5Screate_simple(1, &dim, NULL);
  hid_t dataset= H5Dcreate(h_twist, "twist_angle", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,xk);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  if((xk[0]*xk[0]+xk[1]*xk[1]+xk[2]*xk[2]) < 1e-12)
    is_gamma=1;
  else
    is_gamma=0;
  
  /* create band#/spin# groups so that H5Gopen can be used */
  for(int ib=0; ib<*nband; ib++) 
  {
    char bandname[16];
    sprintf(bandname,"band%i",ib);
    hid_t h_band = H5Gcreate(h_twist,bandname,0);
    for(int ispin=0; ispin<*nspin; ispin++)
    {
      char spinname[16];
      sprintf(spinname,"spin%i",ispin);
      hid_t h_spin = H5Gcreate(h_band,spinname,0);
      H5Gclose(h_spin);
    }
    H5Gclose(h_band);
  }
}
/* close twist# */
void F77_FUNC_(pwhdf_close_twist,PWHDF_CLOSE_TWIST)()
{
  if(h_twist>=0) H5Gclose(h_twist);
  h_twist=-1;
}

/* write eigen value and eigen vector for (ibnd, ispin) */
void F77_FUNC_(pwhdf_write_band,PWHDF_WRITE_BAND)(const int* ibnd,
    const int* ispin, const double* e,
   const double* eigv, const int* ngtot)
{
  char spinname[16];
  sprintf(spinname,"band%i/spin%i",(*ibnd)-1,(*ispin)-1);
  hid_t h_spin = H5Gopen(h_twist,spinname);
  hsize_t dim=1;
  hid_t dataspace= H5Screate_simple(1, &dim, NULL);
  hid_t dataset= H5Dcreate(h_spin, "eigenvalue", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,e);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  hsize_t dims[2];
  dims[0] = *ngtot;
  dims[1] = 2;
  dataspace  = H5Screate_simple(2, dims, NULL);
  dataset =  H5Dcreate(h_spin, "eigenvector", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,eigv);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  H5Gclose(h_spin);
}

/* open mainbody:eigenstates */
void F77_FUNC_(pwhdf_open_eigr,PWHDF_OPEN_EIGR)(const int* nr1, const int* nr2, const int*nr3)
{
  /* swap the order : CHECK THIS! */
  h_ngrid[0]=*nr1;
  h_ngrid[1]=*nr2;
  h_ngrid[2]=*nr3;
  h_ngridtot=(*nr1)*(*nr2)*(*nr3);
  char wfrname[16];
  sprintf(wfrname,"eigenstates_%i_%i_%i",h_ngrid[0],h_ngrid[1],h_ngrid[2]);
  if(h_main>=0) H5Gclose(h_main);
  h_main = H5Gcreate(h_file,wfrname,0);
}

/* close twist# */
void F77_FUNC_(pwhdf_close_eigr,PWHDF_CLOSE_EIGR)()
{
  if(h_main>=0) H5Gclose(h_main);
  h_main=-1;
}
/* write eigen value and eigen vector for (ibnd, ispin) */
void F77_FUNC_(pwhdf_write_wfr,PWHDF_WRITE_WFR)(const int* ibnd,
    const int* ispin, const double* e,
    const double* eigr, const int* use_complex)
{
  char spinname[16];
  sprintf(spinname,"band%i/spin%i",(*ibnd)-1,(*ispin)-1);
  /*sprintf(spinname,"band_%i",(*ibnd)-1);*/
  hid_t h_spin = H5Gopen(h_twist,spinname);

  /* write eigenvalue */
  hsize_t dim=1;
  hid_t dataspace= H5Screate_simple(1, &dim, NULL);
  hid_t dataset= H5Dcreate(h_spin, "eigenvalue", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,e);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  hsize_t dims_out=(*use_complex)?4:3;
  hsize_t dims[4];
  dims[0] = h_ngrid[0];
  dims[1] = h_ngrid[1];
  dims[2] = h_ngrid[2];
  dims[3] = 2;
  dataspace  = H5Screate_simple(dims_out, dims, NULL);
  dataset =  H5Dcreate(h_spin, "eigenvector", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,eigr);

  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Gclose(h_spin);
}

/* write eigen value and eigen vector for (ibnd, ispin) */
void F77_FUNC_(pwhdf_write_rho,PWHDF_WRITE_RHO)(const double* rho, const double* rhog, int ngm)
{
  /* write eigenvector */
  hsize_t dims[3];
  dims[0] = h_ngrid[0];
  dims[1] = h_ngrid[1];
  dims[2] = h_ngrid[2];
  hid_t dataspace  = H5Screate_simple(3, dims, NULL);
  hid_t dataset =  H5Dcreate(h_file, "chargedensity_r", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,rho);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  /*
  hsize_t gdims[2];
  gdims[0]=ngm;
  gdims[1]=2;
  dataspace  = H5Screate_simple(2, gdims, NULL);
  dataset =  H5Dcreate(h_file, "chargedensity_g", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,rhog);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  */

  /* testing with paraview/vtk
  if(is_gamma)
  {
    char vtkname[32];
    sprintf(vtkname,"band%i.vtk",(*ibnd)-1);
    FILE *vtk=fopen(vtkname,"w");

    fprintf(vtk,"# vtk DataFile Version 3.0\n");
    fprintf(vtk,"vtk output\n");
    fprintf(vtk,"ASCII\n");
    fprintf(vtk,"DATASET STRUCTURED_POINTS\n");
    fprintf(vtk,"DIMENSIONS %i %i %i\n",h_ngrid[0],h_ngrid[1],h_ngrid[2]);
    fprintf(vtk,"ORIGIN 0 0 0\n");
    fprintf(vtk,"SPACING 1 1 1\n");
    fprintf(vtk,"\nPOINT_DATA %i\n",h_ngridtot);
    fprintf(vtk,"SCALARS scalars float\n");
    fprintf(vtk,"LOOKUP_TABLE default\n");

    for(int i=0,i2=0; i<h_ngridtot;i+=10)
    { 
      for(int j=0; j<10; j++,i2+=2) fprintf(vtk,"%12.6e ",eigr[i2]*eigr[i2]);
      fprintf(vtk,"\n");
    }
    fprintf(vtk,"\n");
    fclose(vtk);
  }
  */
}
