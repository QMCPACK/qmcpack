//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef OHMMS_HDF_STL_NUMERICATTRIBIO_H
#define OHMMS_HDF_STL_NUMERICATTRIBIO_H
#include "OhmmsData/HDFAttribIO.h"
#include "Numerics/HDFNumericAttrib.h"
#include <complex>
#include <bitset>

namespace qmcplusplus
{
/** Specialization for std::vector<double> */
template<>
struct HDFAttribIO<std::vector<int> >: public HDFAttribIOBase
{

  typedef std::vector<int> ArrayType_t;
  ArrayType_t&  ref;
  bool replace;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a, bool reuse=false):ref(a),replace(reuse) { }

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
      hsize_t dim = ref.size();
      hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
      hid_t dataset =
        H5Dcreate(grp, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
      hid_t ret =
        H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref[0]);
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
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
    hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref[0]);
    H5Sclose(dataspace);
    H5Dclose(h1);
  }

};

template<>
struct HDFAttribIO<std::vector<double> >: public HDFAttribIOBase
{

  typedef std::vector<double> ArrayType_t;
  std::vector<hsize_t> Dim;
  ArrayType_t&  ref;
  bool replace;
  int offset_;

  HDFAttribIO<ArrayType_t>(ArrayType_t& a, bool over=false): ref(a),replace(over),offset_(0)
  {
    Dim.resize(1,a.size());
  }

  HDFAttribIO<ArrayType_t>(ArrayType_t& a, std::vector<int>& dim, int offset=0):
    ref(a), replace(false),offset_(offset)
  {
    Dim.resize(dim.size());
    for(int i=0; i<dim.size(); i++)
      Dim[i]=static_cast<hsize_t>(dim[i]);
  }

  HDFAttribIO<ArrayType_t>(ArrayType_t& a, hsize_t ndim, hsize_t* dim, int offset=0):
    ref(a), replace(false), offset_(offset)
  {
    Dim.resize(ndim);
    for(int i=0; i<ndim; i++)
      Dim[i]= dim[0];
  }

  inline void write(hid_t grp, const char* name)
  {
    //hsize_t dim = ref.size();
    //hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    if(replace)
    {
      hid_t dataset = H5Dopen(grp,name);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(ref[offset_]));
      H5Dclose(dataset);
    }
    else
    {
      hid_t dataspace  = H5Screate_simple(Dim.size(), &Dim[0], NULL);
      hid_t dataset =
        H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      hid_t ret =
        H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref[offset_]);
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

  inline void read(hid_t grp, const char* name)
  {
    hid_t h1 = H5Dopen(grp, name);
    hid_t dataspace = H5Dget_space(h1);
    //hsize_t dims_out[1];
    std::vector<hsize_t> dims_out(Dim);
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n
    = H5Sget_simple_extent_dims(dataspace, &dims_out[0], NULL);
    Dim=dims_out;
    hsize_t ntot=Dim[0];
    for(int i=1; i<Dim.size(); i++)
      ntot*=Dim[i];
    if(ref.size() != int(ntot))
    {
      ref.resize(int(ntot));
    }
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref[0]);
    H5Sclose(dataspace);
    H5Dclose(h1);
  }

};

template<>
struct HDFAttribIO<std::vector<std::complex<double> > >: public HDFAttribIOBase
{
  //NOTE: This specialization assumes each complex number was/will be saved
  // as a pair of doubles. This is checked.

  typedef std::vector<std::complex<double> > ArrayType_t;
  std::vector<hsize_t> Dim;
  ArrayType_t&  ref;

  //Assumes complex stored as a pair of floats/doubles.
  HDFAttribIO<ArrayType_t>(ArrayType_t& a): ref(a)
  {
    Dim.resize(2);
    Dim[0] = a.size();
    Dim[1] = 2;
  }

  inline void write(hid_t grp, const char* name)
  {
    //hsize_t dim = ref.size();
    //hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
    hid_t dataspace  = H5Screate_simple(Dim.size(), &Dim[0], NULL);
    hid_t dataset =
      H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
    hid_t ret =
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref[0]);
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
    hid_t dataspace = H5Dget_space(h1);
    //hsize_t dims_out[1];
    std::vector<hsize_t> dims_out(Dim);
    int rank = H5Sget_simple_extent_ndims(dataspace);
    if(rank!=2)
    {
      //LOGMSG("Error: HDF rank does not match for complex vector")
      return;
    }
    int status_n
    = H5Sget_simple_extent_dims(dataspace, &dims_out[0], NULL);
    if(dims_out[1]!=2)
    {
      //LOGMSG("Error: HDF rank does not match for complex vector")
      return;
    }
    Dim=dims_out;
    hsize_t ntot=Dim[0];
    if(ref.size() != int(ntot))
    {
      ref.resize(int(ntot));
    }
    hid_t ret = H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ref[0]);
    H5Sclose(dataspace);
    H5Dclose(h1);
    H5Eset_auto (func, client_data);
  }

};

template<std::size_t N>
struct HDFAttribIO<std::bitset<N> >: public HDFAttribIOBase
{
  //NOTE: This specialization assumes each complex number was/will be saved
  // as a pair of doubles. This is checked.

  typedef std::bitset<N> ArrayType_t;
  ArrayType_t& ref;
  bool replace;

  //Assumes complex stored as a pair of floats/doubles.
  HDFAttribIO<ArrayType_t>(ArrayType_t& a, bool reuse=false): ref(a), replace(reuse)
  {
  }

  inline void write(hid_t grp, const char* name)
  {
    unsigned long c=ref.to_ulong();
    HDFAttribIO<unsigned long> hc(c,replace);
    hc.write(grp,name);
  }

  inline void read(hid_t grp, const char* name)
  {
    unsigned long c=ref.to_ulong();
    HDFAttribIO<unsigned long> hc(c);
    hc.read(grp,name);
    ref=c;
  }
};
}
#endif
