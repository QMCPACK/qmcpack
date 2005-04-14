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
#ifndef OHMMS_HDF_NUMERICATTRIBIO_H
#define OHMMS_HDF_NUMERICATTRIBIO_H
#include "OhmmsData/HDFAttribIO.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#ifdef HAVE_LIBBLITZ
#include <blitz/array.h>
#endif

/** Specialization for hsize_t */
template<>
struct HDFAttribIO<hsize_t>: public HDFAttribIOBase {

  hsize_t& ref;

  HDFAttribIO<hsize_t>(hsize_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name) {
    hsize_t dim = 1;
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  
      H5Dcreate(grp, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref);
    H5Sclose(dataspace);

    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name) {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref);
    H5Dclose(h1);
  }
};

/*

template<unsigned D>
struct HDFAttribIO<TinyVector<double,D> >: public HDFAttribIOBase {

typedef TinyVector<double,D> data_type;

data_type& ref;

HDFAttribIO<data_type>(data_type& a):ref(a) { }

inline void write(hid_t grp, const char* name) {
hsize_t dim = D;
hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
hid_t dataset =  
H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
hid_t ret = 
H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0]));
H5Sclose(dataspace);
H5Dclose(dataset);
}

inline void read(hid_t  grp, const char* name) {
hid_t h1 = H5Dopen(grp, name);
hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(ref[0]));
H5Dclose(h1);
}
};

*/

/** Specialization for Vector<double> */
template<>
struct HDFAttribIO<Vector<double> >: public HDFAttribIOBase {

  typedef Vector<double> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name) {

    hsize_t dim = ref.size();
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
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
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Sclose(dataspace);
    H5Dclose(h1);
  }

};

/** Specialization for Vector<int>  */
template<>
struct HDFAttribIO<Vector<int> >: public HDFAttribIOBase {

  typedef Vector<int> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name) {

    hsize_t dim = ref.size();
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  
      H5Dcreate(grp, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
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
    hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Sclose(dataspace);
    H5Dclose(h1);
  }

};

/** Specialization for Matrix<double> */
template<>
struct HDFAttribIO<Matrix<double> >: public HDFAttribIOBase {

  typedef Matrix<double> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name) {
    hsize_t dim = ref.size();
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name) {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Dclose(h1);
  }

};

#ifdef HAVE_LIBBLITZ
/** Specialization for blitz::Array<TinyVector<double,D>,2> */
template<unsigned D>
struct HDFAttribIO<blitz::Array<TinyVector<double,D>,2> >: public HDFAttribIOBase {

  typedef blitz::Array<TinyVector<double,D>,2> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }
 
  inline void write(hid_t grp, const char* name) {
    int rank = 3;
    hsize_t dim[rank];
    dim[0] = ref.extent(0);
    dim[1] = ref.extent(1);
    dim[2] = D;
    hid_t dataspace  = H5Screate_simple(rank, dim, NULL);
    hid_t dataset =  
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name) {
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[3];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if((ref.extent(0) != (unsigned long)dims_out[0]) || (ref.extent(1) != (unsigned long)dims_out[1])){
      //   cout << "dimensions not equal" << endl;
      ref.resize(dims_out[0],dims_out[1]);
    }
       
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Dclose(h1);
  }
  
};
#endif

/** Specialization for Vector<TinyVector<double,D> > */
template<unsigned D>
struct HDFAttribIO<Vector<TinyVector<double,D> > >: public HDFAttribIOBase {

  typedef Vector<TinyVector<double,D> > ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }
  
  inline void write(hid_t grp, const char* name) {
    hsize_t dim[2] = {ref.size(), D};
    hid_t dataspace  = H5Screate_simple(2, dim, NULL);
    hid_t dataset =  
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret = 
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name) {
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[2];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if(ref.size() != (unsigned long)dims_out[0]){
      //   cout << "dimensions not equal" << endl;
      ref.resize(dims_out[0]);
    }
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Dclose(h1);
  }

};

/** Specialization for string */
template<>
struct HDFAttribIO<string>: public HDFAttribIOBase {

  typedef string ArrayType_t;
  ArrayType_t& ref;
  hid_t str80; 

  HDFAttribIO<ArrayType_t>(ArrayType_t& a): ref(a) {
    str80 = H5Tcopy(H5T_C_S1);
    H5Tset_size(str80,80);
  }
  
  inline void write(hid_t grp, const char* name) {

    hsize_t dim = 1;
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =  
      //    H5Dcreate(grp, name, H5T_NATIVE_CHAR, dataspace, H5P_DEFAULT);
      H5Dcreate(grp, name, str80, dataspace, H5P_DEFAULT);
    hid_t ret = 
      //   H5Dwrite(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
      H5Dwrite(dataset, str80, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name) {
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[1];
    //    int rank = H5Sget_simple_extent_ndims(dataspace);
    //    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    // if(ref.size() != (unsigned long)dims_out[0]){
    //   cout << "dimensions not equal" << endl;
    //    ref.resize(dims_out[0]);
    //  }
    hid_t ret = H5Dread(h1, str80, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref);
    H5Dclose(h1);
  }

};



#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
