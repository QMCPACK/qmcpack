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
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_HDF_TRAITS_H
#define QMCPLUSPLUS_HDF_TRAITS_H
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

/** common traits for atomic data type
 *
 * Apply to intrinsic types, e.g, double, int
 */
template<typename T>
struct atomic_type_traits
{
  enum {is_array_type=0};
  enum {array_dimension=1};
  typedef T* pointer_type;
  typedef const T* const_pointer_type;
  hsize_t dims[1];
  hid_t h5_data_type;
  T& ref_;
  atomic_type_traits(T& a):ref_(a)
  {
    dims[0]=1;
  }
  inline hsize_t size() const
  {
    return 1;
  }
  inline const_pointer_type data() const
  {
    return &ref_;
  }
  inline pointer_type data()
  {
    return &ref_;
  }
};

/** dummy class declaration  for h5 datatype traits
 *
 * To be specialized. h5_data_type has to be set correctly for T
 * In order to be used for HDF, each specialization should provide
 * - array_dimension
 * - dims[array_dimension or larger]
 * - h5_data_type : H5T_NATIVE_* for the data type
 */
template<typename T> struct h5_data_type_traits
  { };

/** specialization for double */
template<>
struct  h5_data_type_traits<double>: public atomic_type_traits<double>
{
  h5_data_type_traits<double>(double& a):atomic_type_traits<double>(a)
  {
    h5_data_type=H5T_NATIVE_DOUBLE;
  }
};

/** specialization for int */
template<>
struct h5_data_type_traits<int>: public atomic_type_traits<int>
{
  h5_data_type_traits<int>(int& a):atomic_type_traits<int>(a)
  {
    h5_data_type=H5T_NATIVE_INT;
  }
};

/** free function to return the H5T_NATIVE* for T
 */
template<typename T>
inline hid_t get_h5_data_type(T x)
{
  h5_data_type_traits<T> dummy(x);
  return dummy.h5_data_type;
}


/* specialization for vector<T>
 */
template<typename T>
struct h5_data_type_traits<vector<T> >
{
  enum {is_array_type=1};
  enum {array_dimension=1};
  hsize_t dims[1];
  hid_t h5_data_type;
  vector<T>& ref_;
  h5_data_type_traits(vector<T>& a):ref_(a)
  {
    dims[0]=ref_.size();
    h5_data_type=get_h5_data_type(T());
  }
  inline hsize_t size() const
  {
    return ref_.size();
  }
  inline const T* data() const
  {
    return &ref_[0];
  }
  inline T* data()
  {
    return &ref_[0];
  }
};

/* specialization for TinyVector<T,D>
 */
template<typename T, unsigned D>
struct h5_data_type_traits<TinyVector<T,D> >
{
  enum {is_array_type=1};
  enum {array_dimension=1};
  hsize_t dims[1];
  hid_t h5_data_type;
  TinyVector<T,D>& ref_;

  h5_data_type_traits(TinyVector<T,D>& a):ref_(a)
  {
    dims[0]=D;
    h5_data_type=get_h5_data_type(T());
  }
  inline hsize_t size() const
  {
    return D;
  }
  inline const T* data() const
  {
    return &ref_[0];
  }
  inline T* data()
  {
    return &ref_[0];
  }
};

/* specialization for vector<TinyVector<T,D> >
 */
template<typename T, unsigned D>
struct h5_data_type_traits<vector<TinyVector<T,D> > >
{
  enum {is_array_type=1};
  enum {array_dimension=2};
  hsize_t dims[3];
  hid_t h5_data_type;
  vector<TinyVector<T,D> >& ref_;

  h5_data_type_traits(vector<TinyVector<T,D> >& a):ref_(a)
  {
    dims[0]=ref_.size();
    dims[1]=D;
    dims[2]=D*ref_.size();
    h5_data_type=get_data_type(T());
  }

  inline hsize_t size() const
  {
    return dims[2];
  }
  inline const T* data() const
  {
    return &ref_[0];
  }
  inline T* data()
  {
    return &ref_[0];
  }
};


template<typename T, unsigned D>
struct h5_data_type_traits<Array<T,D> >
{
  enum {is_array_type=1};
  enum {array_dimension=D};
  hsize_t dims[D];
  hid_t h5_data_type;
  Array<T,D>& ref_;

  h5_data_type_traits(Array<T,D>& a):ref_(a)
  {
    for(int i=0; i<D; ++i)
    {
      dims[i]=ref_.size(i);
    }
    h5_data_type=get_data_type(T());
  }
  inline hsize_t size() const
  {
    hsize_t tot=dims[1];
    for(int i=1; i<D; ++i)
      tot*=dims[i];
    return tot;
  }
  inline const T* data() const
  {
    return ref_->data();
  }
  inline T* data()
  {
    return ref_->data();
  }
};

template<typename T>
struct HDFSource: public h5_data_type_traits<T>//, public HDFAttribIOBase
{
  using h5_data_type_traits<T>::is_array_type;
  using h5_data_type_traits<T>::array_dimension;
  using h5_data_type_traits<T>::h5_data_type;

  bool overwrite;
  HDFSource(T& a, bool replace=false):h5_data_type_traits<T>(a),overwrite(replace) {}

  inline void write(hid_t grp, const char* name)
  {
    hid_t dataspace  = H5Screate_simple(array_dimension, dims, NULL);
    hid_t dataset =  H5Dcreate(grp, name, datatype, dataspace, H5P_DEFAULT);
    hid_t ret = H5Dwrite(dataset, h5_data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ref);
    H5Sclose(dataspace);
  }

  inline void read()
  {
    //cout << "h5 datatype = " << this->h5_data_type() << " " << this->ref_ << endl;
    cout << "h5 datatype = " << h5_data_type << endl;
    cout << "address = " << this->data() << endl;
    if(is_array_type)
    {
      cout << "This is an array. Its size = " << this->size() << endl;
    }
    else
    {
      cout << "This is not an array " << endl;
    }
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2764 $   $Date: 2008-06-26 10:21:31 -0500 (Thu, 26 Jun 2008) $
 * $Id: HDFNumericAttrib.h 2764 2008-06-26 15:21:31Z jnkim $
 ***************************************************************************/
