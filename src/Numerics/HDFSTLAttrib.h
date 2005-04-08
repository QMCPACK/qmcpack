//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
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
#ifndef OHMMS_HDF_STL_NUMERICATTRIBIO_H
#define OHMMS_HDF_STL_NUMERICATTRIBIO_H
#include "OhmmsData/HDFAttribIO.h"

/** Specialization for std::vector<double> */
template<>
struct HDFAttribIO<std::vector<double> >: public HDFAttribIOBase {

  typedef std::vector<double> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name) {

    hsize_t dim = ref.size();
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref[0]);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name) {

    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[1];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if(ref.size() != int(dims_out[0])){
      ref.resize(int(dims_out[0]));
    }
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref[0]);
    H5Sclose(dataspace);
    H5Dclose(h1);
  }

};

template<>
struct HDFAttribIO<std::vector<int> >: public HDFAttribIOBase {

  typedef std::vector<int> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name) {

    hsize_t dim = ref.size();
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  
      H5Dcreate(grp, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref[0]);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name) {

    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[1];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if(ref.size() != int(dims_out[0])){
      ref.resize(int(dims_out[0]));
    }
    hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref[0]);
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
