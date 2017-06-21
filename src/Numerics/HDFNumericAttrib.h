//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef OHMMS_HDF_NUMERICATTRIBIO_H
#define OHMMS_HDF_NUMERICATTRIBIO_H
#include "OhmmsData/HDFAttribIO.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#ifdef HAVE_LIBBLITZ
#include <blitz/array.h>
#else
#include "OhmmsPETE/OhmmsArray.h"
#endif


namespace qmcplusplus
{


/** Specialization for std::string */
/*
template<>
struct HDFAttribIO<std::string>: public HDFAttribIOBase {

  std::string& ref;

  HDFAttribIO<std::string>( std::string& a):ref(a) { }

  inline void write(hid_t grp, const char* name) {
    hsize_t dim = 1;
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t typeID = H5Tcopy (H5T_C_S1);
    H5Tset_size(typeID, ref.size());
    hid_t dataset =
      H5Dcreate(grp, name, typeID, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, typeID, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0]));
    H5Sclose(dataspace);
    H5Dclose(dataset);
    H5Tclose(typeID);
  }

  inline void read(hid_t grp, const char* name) {
    hid_t h1 = H5Dopen(grp, name);
    hid_t typeID = H5Dget_type(h1);
    hsize_t len = H5Tget_size(typeID);
    ref.resize(len);
    hid_t ret = H5Dread(h1, typeID, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(ref[0]));
    H5Tclose(typeID);
    H5Dclose(h1);
  }
  };*/


/** Specialization for hsize_t */
template<>
struct HDFAttribIO<hsize_t>: public HDFAttribIOBase
{

  hsize_t& ref;

  HDFAttribIO<hsize_t>(hsize_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
    hsize_t dim = 1;
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref);
    H5Dclose(h1);
  }
};

template<>
struct HDFAttribIO<unsigned long>: public HDFAttribIOBase
{

  unsigned long& ref;
  bool replace;

  HDFAttribIO<unsigned long>(unsigned long& a, bool reuse=false):ref(a),replace(reuse) { }

  inline void write(hid_t grp, const char* name)
  {
    if(replace)
    {
      hid_t h1 = H5Dopen(grp, name);
      hid_t ret = H5Dwrite(h1, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref);
      H5Dclose(h1);
    }
    else
    {
      hsize_t dim = 1;
      hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
      hid_t dataset =  H5Dcreate(grp, name, H5T_NATIVE_ULONG, dataspace, H5P_DEFAULT);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref);
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref);
    H5Dclose(h1);
  }
};

/** Specialization for int */
template<>
struct HDFAttribIO<int>: public HDFAttribIOBase
{
  int& ref;
  bool replace;

  HDFAttribIO<int>(int& a, bool reuse=false):ref(a),replace(reuse) { }

  inline void write(hid_t grp, const char* name)
  {
    if(replace)
      //herr_t status = H5Eset_auto(NULL, NULL);
      //status = H5Gget_objinfo (grp, name, 0, NULL);
      //if(status == 0)
    {
      hid_t h1 = H5Dopen(grp, name);
      hid_t ret = H5Dwrite(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref);
      H5Dclose(h1);
    }
    else
    {
      hsize_t dim = 1;
      hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
      hid_t dataset =  H5Dcreate(grp, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref);
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

  inline void read(hid_t grp, const char* name)
  {
    // Turn off error printing
    H5E_auto_t func;
    void *client_data;
    H5Eget_auto (&func, &client_data);
    H5Eset_auto (NULL, NULL);
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ref);
    H5Dclose(h1);
    H5Eset_auto (func, client_data);
  }

};


/** Specialization for double */
template<>
struct HDFAttribIO<double>: public HDFAttribIOBase
{

  double& ref;

  HDFAttribIO<double>(double& a):ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
    hsize_t dim = 1;
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name)
  {
    H5E_auto_t func;
    void *client_data;
    H5Eget_auto (&func, &client_data);
    H5Eset_auto (NULL, NULL);
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = -1;
    if (h1 >= 0)
    {
      ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref);
      H5Dclose(h1);
    }
    H5Eset_auto (func, client_data);
    // return (ret >= 0);
  }
};

template<unsigned D>
struct HDFAttribIO<TinyVector<double,D> >: public HDFAttribIOBase
{

  typedef TinyVector<double,D> data_type;

  data_type& ref;
  bool replace;

  inline HDFAttribIO<data_type>(data_type& a, bool over=false):ref(a),replace(over) { }

  inline void write(hid_t grp, const char* name)
  {
    if(replace)
    {
      hid_t dataset = H5Dopen(grp,name);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0]));
      H5Dclose(dataset);
    }
    else
    {
      hsize_t dim = D;
      hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
      hid_t dataset =  H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0]));
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

  inline void read(hid_t  grp, const char* name)
  {
    H5E_auto_t func;
    void *client_data;
    H5Eget_auto (&func, &client_data);
    H5Eset_auto (NULL, NULL);
    hid_t h1 = H5Dopen(grp, name);
    if(h1>-1)
    {
      hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(ref[0]));
      H5Dclose(h1);
    }
    H5Eset_auto (func, client_data);
  }
};

template<unsigned D>
struct HDFAttribIO<Tensor<double,D> >: public HDFAttribIOBase
{

  typedef Tensor<double,D> data_type;

  data_type& ref;
  bool replace;

  inline HDFAttribIO<data_type>(data_type& a, bool over=false):ref(a),replace(over) { }

  inline void write(hid_t grp, const char* name)
  {
    if(replace)
    {
      hid_t dataset = H5Dopen(grp,name);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref(0,0)));
      H5Dclose(dataset);
    }
    else
    {
      hsize_t dim[2] = {D, D};
      hid_t dataspace  = H5Screate_simple(2, dim, NULL);
      hid_t dataset =  H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref(0,0)));
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

  inline void read(hid_t  grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(ref(0,0)));
    H5Dclose(h1);
  }
};



template<unsigned D>
struct HDFAttribIO<TinyVector<int,D> >: public HDFAttribIOBase
{

  typedef TinyVector<int,D> data_type;
  data_type& ref;
  bool replace;

  HDFAttribIO<data_type>(data_type& a, bool reuse=false):ref(a),replace(reuse) { }

  inline void write(hid_t grp, const char* name)
  {
    if(replace)
    {
      hid_t dataset = H5Dopen(grp,name);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0]));
      H5Dclose(dataset);
    }
    else
    {
      hsize_t dim = D;
      hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
      hid_t dataset =
        H5Dcreate(grp, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
      hid_t ret =
        H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0]));
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

  inline void read(hid_t  grp, const char* name)
  {
    // Turn off error printing
    H5E_auto_t func;
    void *client_data;
    H5Eget_auto (&func, &client_data);
    H5Eset_auto (NULL, NULL);
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(ref[0]));
    H5Dclose(h1);
    H5Eset_auto (func, client_data);
  }
};

//Specialization for std::vector<TinyVector<int,D> >
template<unsigned D>
struct HDFAttribIO<std::vector<TinyVector<int,D> > >: public HDFAttribIOBase
{

  typedef std::vector<TinyVector<int,D> > data_type;

  data_type& ref;

  HDFAttribIO<data_type>(data_type& a):ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
    hsize_t dim[2] = {ref.size(), D};
    hid_t dataspace  = H5Screate_simple(2, dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0][0]));
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t  grp, const char* name)
  {
    H5E_auto_t func;
    void *client_data;
    H5Eget_auto (&func, &client_data);
    H5Eset_auto (NULL, NULL);
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = -1;
    if (h1 >= 0)
    {
      hid_t dataspace = H5Dget_space(h1);
      hsize_t dim[2];
      int rank = H5Sget_simple_extent_ndims(dataspace);
      int status_n = H5Sget_simple_extent_dims(dataspace, dim, NULL);
      //Resize storage if not equal
      if(ref.size() != (unsigned long)dim[0])
      {
        ref.resize(dim[0]);
      }
      hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(ref[0][0]));
      H5Dclose(h1);
    }
    H5Eset_auto (func, client_data);
  }
};

//Specialization for std::vector<TinyVector<double,D> >
template<unsigned D>
struct HDFAttribIO<std::vector<TinyVector<double,D> > >: public HDFAttribIOBase
{

  typedef std::vector<TinyVector<double,D> > data_type;

  data_type& ref;

  HDFAttribIO<data_type>(data_type& a):ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
    hsize_t dim[2] = {ref.size(), D};
    hid_t dataspace  = H5Screate_simple(2, dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0][0]));
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t  grp, const char* name)
  {
    // Turn off error printing
    H5E_auto_t func;
    void *client_data;
    H5Eget_auto (&func, &client_data);
    H5Eset_auto (NULL, NULL);
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = -1;
    if (h1 >= 0)
    {
      hid_t dataspace = H5Dget_space(h1);
      hsize_t dim[2]= {0,D};
      int rank = H5Sget_simple_extent_ndims(dataspace);
      int status_n = H5Sget_simple_extent_dims(dataspace, dim, NULL);
      //Resize storage if not equal
      if(ref.size() != (unsigned long)dim[0])
      {
        ref.resize(dim[0]);
      }
      ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[0][0]));
      H5Dclose(h1);
    }
    H5Eset_auto (func, client_data);
  }
};


/** Specialization for Vector<double> */
template<>
struct HDFAttribIO<Vector<double> >: public HDFAttribIOBase
{

  typedef Vector<double> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
    hsize_t dim = ref.size();
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[1];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if(ref.size() != int(dims_out[0]))
    {
      ref.resize(int(dims_out[0]));
    }
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Sclose(dataspace);
    H5Dclose(h1);
  }

};


template<>
struct HDFAttribIO<Vector<std::complex<double> > >: public HDFAttribIOBase
{

  typedef Vector<std::complex<double> > ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
    hsize_t dim[2];
    dim[0] = ref.size();
    dim[1] = 2;
    hid_t dataspace  = H5Screate_simple(2, dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[2];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank != 2)
    {
      fprintf (stderr, "Error reading Vector<complex> in HDFNumericAttrib.h.\n");
      abort();
    }
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if(ref.size() != int(dims_out[0]))
    {
      ref.resize(int(dims_out[0]));
    }
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Sclose(dataspace);
    H5Dclose(h1);
  }

};



/** Specialization for Vector<int>  */
template<>
struct HDFAttribIO<Vector<int> >: public HDFAttribIOBase
{

  typedef Vector<int> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
    hsize_t dim = ref.size();
    hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[1];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if(ref.size() != int(dims_out[0]))
    {
      ref.resize(int(dims_out[0]));
    }
    hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Sclose(dataspace);
    H5Dclose(h1);
  }

};

/** Specialization for Vector<TinyVector<double,D> > */
template<unsigned D>
struct HDFAttribIO<Vector<TinyVector<double,D> > >: public HDFAttribIOBase
{

  typedef Vector<TinyVector<double,D> > ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
    hsize_t dim[2] = {ref.size(), D};
    hid_t dataspace  = H5Screate_simple(2, dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[2];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if(ref.size() != (unsigned long)dims_out[0])
    {
      //   std::cout << "dimensions not equal" << std::endl;
      ref.resize(dims_out[0]);
    }
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Dclose(h1);
  }

};


/** Specialization for Matrix<int> */
template<>
struct HDFAttribIO<Matrix<int> >: public HDFAttribIOBase
{

  typedef Matrix<int> ArrayType_t;
  ArrayType_t&  ref;
  bool replace;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a, bool reuse=false):ref(a),replace(reuse) { }

  inline void write(hid_t grp, const char* name)
  {
    if(replace)
      //  herr_t status = H5Eset_auto(NULL, NULL);
      //  status = H5Gget_objinfo (grp, name, 0, NULL);
      //  if(status == 0)
    {
      hid_t h1 = H5Dopen(grp, name);
      hid_t ret= H5Dwrite(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
      H5Dclose(h1);
    }
    else
    {
      hsize_t dims[2];
      dims[0]=ref.rows();
      dims[1]=ref.cols();
      hid_t dataspace  = H5Screate_simple(2, dims, NULL);
      hid_t dataset =
        H5Dcreate(grp, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
      hid_t ret =
        H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Dclose(h1);
  }
};


/** Specialization for Matrix<double> */
template<>
struct HDFAttribIO<Matrix<double> >: public HDFAttribIOBase
{

  typedef Matrix<double> ArrayType_t;
  ArrayType_t&  ref;
  bool replace;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a, bool reuse=false):ref(a),replace(reuse) { }

  inline void write(hid_t grp, const char* name)
  {
    if(replace)
      //  herr_t status = H5Eset_auto(NULL, NULL);
      //  status = H5Gget_objinfo (grp, name, 0, NULL);
      //  if(status == 0)
    {
      hid_t h1 = H5Dopen(grp, name);
      hid_t ret= H5Dwrite(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
      H5Dclose(h1);
    }
    else
    {
      hsize_t dims[2];
      dims[0]=ref.rows();
      dims[1]=ref.cols();
      hid_t dataspace  = H5Screate_simple(2, dims, NULL);
      hid_t dataset =
        H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      hid_t ret =
        H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Dclose(h1);
  }
};

/** Specialization for Matrix<TinyVector<double,D>> */
template<unsigned D>
struct HDFAttribIO<Matrix<TinyVector<double,D> > >: public HDFAttribIOBase
{

  typedef TinyVector<double,D> Component_t;
  typedef Matrix<Component_t>  ArrayType_t;
  ArrayType_t&  ref;
  bool replace;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a, bool reuse=false):ref(a), replace(reuse) { }

  inline void write(hid_t grp, const char* name)
  {
    const int rank = 3;
    hsize_t dims[rank], maxdims[rank];
    dims[0] = ref.rows();
    dims[1] = ref.cols();
    dims[2] = D;
    maxdims[0]=H5S_UNLIMITED;
    maxdims[1] = ref.cols();
    maxdims[2] = D;
    if(replace)
    {
      hid_t dataset =  H5Dopen(grp, name);
      herr_t status=H5Dextend(dataset,dims);
      hid_t memspace = H5Screate_simple(rank, dims, NULL);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, H5S_ALL, H5P_DEFAULT,ref.data());
      H5Sclose(memspace);
      H5Dclose(dataset);
    }
    else
    {
      hid_t dataspace  = H5Screate_simple(rank, dims, maxdims);
      hid_t p = H5Pcreate (H5P_DATASET_CREATE);
      H5Pset_chunk(p,rank,dims);
      hid_t dataset =  H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, p);
      hid_t memspace = H5Screate_simple(rank, dims, NULL);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,ref.data());
      H5Sclose(memspace);
      H5Sclose(dataspace);
      H5Dclose(dataset);
      H5Pclose(p);
    }
    /*
    if(replace)
    {
      hid_t dataset =  H5Dopen(grp, name);
      hid_t ret =
        H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
      H5Dclose(dataset);
    }
    else
    {
      const int rank = 3;
      hsize_t dims[rank];
      dims[0] = ref.rows(); dims[1] = ref.cols(); dims[2] = D;
      hid_t dataspace  = H5Screate_simple(rank, dims, NULL);
      hid_t dataset =
        H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      hid_t ret =
        H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
    */
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[3];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if((ref.rows() != int(dims_out[0])) || (ref.cols() != int(dims_out[1])))
    {
      ref.resize(dims_out[0],dims_out[1]);
    }
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Dclose(h1);
  }
};

#ifdef HAVE_LIBBLITZ
/** Specialization for blitz::Array<TinyVector<double,D>,2> */
template<unsigned D>
struct HDFAttribIO<blitz::Array<TinyVector<double,D>,2> >: public HDFAttribIOBase
{

  typedef blitz::Array<TinyVector<double,D>,2> ArrayType_t;
  ArrayType_t&  ref;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }

  inline void write(hid_t grp, const char* name)
  {
    const int rank = 3;
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

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[3];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    if((ref.extent(0) != (unsigned long)dims_out[0]) || (ref.extent(1) != (unsigned long)dims_out[1]))
    {
      //   std::cout << "dimensions not equal" << std::endl;
      ref.resize(dims_out[0],dims_out[1]);
    }
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
    H5Dclose(h1);
  }

};
#else
/** Specialization for Array<double,D> */
template <unsigned D>
struct HDFAttribIO<Array<double,D> >: public HDFAttribIOBase
{
  typedef Array<double,D> ArrayType_t;
  ArrayType_t&  ref;
  HDFAttribIO<ArrayType_t>(ArrayType_t& a):ref(a) { }
  inline void write(hid_t grp, const char* name)
  {
    const int rank = D;
    hsize_t dim[rank];
    for(int i=0; i<rank; ++i)
      dim[i]=ref.size(0);
    hid_t dataspace  = H5Screate_simple(rank, dim, NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }

  inline void read(hid_t grp, const char* name)
  {
    // Turn off error printing
    H5E_auto_t func;
    void *client_data;
    H5Eget_auto (&func, &client_data);
    H5Eset_auto (NULL, NULL);
    hid_t h1 = H5Dopen(grp, name);
    if (h1 > 0)
    {
      hid_t dataspace = H5Dget_space(h1);
      hsize_t dims[D];
      H5Sget_simple_extent_dims(dataspace, dims, NULL);
      ref.resize(dims);//[0], dims[1], dims[2]);
      hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
      H5Sclose(dataspace);
      H5Dclose(h1);
    }
    H5Eset_auto (func, client_data);
  }
};

/** Specialization for Array<std::complex<double>,D> */
template<unsigned D>
struct HDFAttribIO<Array<std::complex<double>,D> >: public HDFAttribIOBase
{
  typedef Array<std::complex<double>,D> ArrayType_t;

  ArrayType_t&  ref;
  bool replace;
  HDFAttribIO<ArrayType_t>(ArrayType_t& a, bool reuse=false):ref(a), replace(reuse) { }
  inline void write(hid_t grp, const char* name)
  {
    const int rank = D+1;
    hsize_t dim[rank];
    for(int i=0; i<D; ++i)
      dim[i]=ref.size(i);
    dim[D] = 2;
    if(replace)
    {
      hid_t dataset =  H5Dopen(grp, name);
      herr_t status = H5Dextend(dataset,dim);
      hid_t dataspace = H5Screate_simple(rank, dim, NULL);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, H5S_ALL, H5P_DEFAULT,ref.data());
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
    else
    {
      hid_t dataspace  = H5Screate_simple(rank, dim, NULL);
      hid_t dataset =
        H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      hid_t ret =
        H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,ref.data());
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    if(h1>=-1)
    {
      hid_t dataspace = H5Dget_space(h1);
      hsize_t dims[D+1];
      H5Sget_simple_extent_dims(dataspace, dims, NULL);
      ref.resize(dims);
      hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref.data());
      H5Sclose(dataspace);
      H5Dclose(h1);
    }
  }
};
#endif
}
#endif
