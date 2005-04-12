//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_HDF_TRICUBICSPLINE_H
#define OHMMS_HDF_TRICUBICSPLINE_H

#include "OhmmsData/HDFAttribIO.h"
#include "Numerics/TriCubicSplineT.h"

/** Specialization for std::vector<double> */
template<>
struct HDFAttribIO<TriCubicSplineT<double> >: public HDFAttribIOBase {

  typedef TriCubicSplineT<double> Data_t;
  Data_t&  ref;

  HDFAttribIO<Data_t>(Data_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name) {
    hsize_t dim[4];
    dim[0]=ref.nX; dim[1]=ref.nY; dim[2]=ref.nZ; dim[3]=8; 
    hid_t dataspace  = H5Screate_simple(4, dim, NULL);
    hid_t dataset =  
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name) {
    hsize_t dim[4];
    dim[0]=ref.nX; dim[1]=ref.nY; dim[2]=ref.nZ; dim[3]=8; 
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Sclose(dataspace);
    H5Dclose(h1);
  }

};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
