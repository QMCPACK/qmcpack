//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF_DATASPACE_TRAITS_H
#define QMCPLUSPLUS_HDF_DATASPACE_TRAITS_H
/**@file hdf_dataspace.h
 * @brief define h5_space_type to handle basic datatype for hdf5
 *
 * h5_space_type is a helper class for h5data_proxy and used internally
 * - h5_space_type<T,DS>
 * - h5_space_type<std::complex<T>,DS>
 * - h5_space_type<TinyVector<T,D>,DS>
 * - h5_space_type<TinyVector<std::complex<T>,D>,DS>
 * - h5_space_type<Tensor<T,D>,DS>
 * - h5_space_type<Tensor<std::complex<T>,D>,DS>
 */


#include <io/hdf_datatype.h>
#include <complex>
#include <OhmmsPETE/TinyVector.h>
#include <OhmmsPETE/Tensor.h>

namespace qmcplusplus
{
/** default struct to define a h5 dataspace, any intrinsic type T
 *
 * \tparm T intrinsic datatype
 * \tparm D dimension of the h5dataspace
 */
template<typename T, hsize_t D>
struct h5_space_type
{
  ///pointer type
  typedef T* pointer;
  ///shape of the dataspace
  hsize_t dims[D];
  ///size, dimension,  of the dataspace
  inline int size() const
  {
    return D;
  }
  ///return the address
  inline pointer get_address(T* a)
  {
    return a;
  }
};

/** specialization of h5_space_type for std::complex<T>
 *
 * Raize the dimension of the space by 1 and set the last dimension=2
 */
template<typename T, hsize_t D>
struct h5_space_type<std::complex<T>,D>
{
  typedef T* pointer;
  hsize_t dims[D+1];
  inline h5_space_type()
  {
    dims[D]=2;
  }
  inline int size() const
  {
    return D+1;
  }
  inline pointer get_address(std::complex<T>* a)
  {
    return reinterpret_cast<T*>(a);
  }
};

/** specialization of h5_space_type for TinyVector<T,D> for any intrinsic type T
 */
template<typename T, unsigned D, hsize_t DS>
struct h5_space_type<TinyVector<T,D>, DS>
{
  typedef T* pointer;
  hsize_t dims[DS+1];
  inline h5_space_type()
  {
    dims[DS]=D;
  }
  inline int size() const
  {
    return DS+1;
  }
  inline pointer get_address(TinyVector<T,D>* a)
  {
    return a->data();
  }
};

/** specialization of h5_space_type for TinyVector<std::complex<T>,D> for std::complex<T>
 */
template<typename T, unsigned D, hsize_t DS>
struct h5_space_type<TinyVector<std::complex<T>,D>, DS>
{
  typedef T* pointer;
  hsize_t dims[DS+2];
  inline h5_space_type()
  {
    dims[DS]=D;
    dims[DS+1]=2;
  }
  inline int size() const
  {
    return DS+2;
  }
  inline pointer get_address(TinyVector<std::complex<T>,D>* a)
  {
    return reinterpret_cast<T*>(a->data());
  }
};

/** specialization of h5_space_type for TinyVector<T,D> for any intrinsic type T
 */
template<typename T, unsigned D, hsize_t DS>
struct h5_space_type<Tensor<T,D>, DS>
{
  typedef T* pointer;
  hsize_t dims[DS+2];
  inline h5_space_type()
  {
    dims[DS]=D;
    dims[DS+1]=D;
  }
  inline int size() const
  {
    return DS+2;
  }
  inline pointer get_address(Tensor<T,D>* a)
  {
    return a->data();
  }
};

/** specialization of h5_space_type for TinyVector<std::complex<T>,D> for std::complex<T>
 */
template<typename T, unsigned D, hsize_t DS>
struct h5_space_type<Tensor<std::complex<T>,D>, DS>
{
  typedef T* pointer;
  hsize_t dims[DS+3];
  inline h5_space_type()
  {
    dims[DS]=D;
    dims[DS+1]=D;
    dims[DS+2]=2;
  }
  inline int size() const
  {
    return DS+2;
  }
  inline pointer get_address(Tensor<std::complex<T>,D>* a)
  {
    return reinterpret_cast<T*>(a->data());
  }
};
}
#endif
