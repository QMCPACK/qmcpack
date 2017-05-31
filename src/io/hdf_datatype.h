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


#ifndef QMCPLUSPLUS_H5DATATYPE_DEFINE_H
#define QMCPLUSPLUS_H5DATATYPE_DEFINE_H

#include <type_traits/scalar_traits.h>
#if defined(HAVE_LIBHDF5)
#include <hdf5.h>
#endif

namespace qmcplusplus
{
#if defined(HAVE_LIBHDF5)
template <typename T>
inline hid_t
get_h5_datatype(const T&)
{
  return H5T_NATIVE_CHAR;
}

#define BOOSTSUB_H5_DATATYPE(CppType, H5DTYPE)                         \
template<>                                                             \
inline hid_t                                                           \
get_h5_datatype< CppType >(const CppType&) { return H5DTYPE; }

BOOSTSUB_H5_DATATYPE(short, H5T_NATIVE_SHORT);
BOOSTSUB_H5_DATATYPE(int, H5T_NATIVE_INT);
BOOSTSUB_H5_DATATYPE(long, H5T_NATIVE_LONG);
BOOSTSUB_H5_DATATYPE(unsigned char, H5T_NATIVE_UCHAR);
BOOSTSUB_H5_DATATYPE(unsigned int, H5T_NATIVE_UINT);
BOOSTSUB_H5_DATATYPE(unsigned long, H5T_NATIVE_ULONG);
//BOOSTSUB_H5_DATATYPE(uint32_t, H5T_NATIVE_UINT32);
//BOOSTSUB_H5_DATATYPE(uint64_t, H5T_NATIVE_UINT64);
BOOSTSUB_H5_DATATYPE(float, H5T_NATIVE_FLOAT);
BOOSTSUB_H5_DATATYPE(double, H5T_NATIVE_DOUBLE);
BOOSTSUB_H5_DATATYPE(std::complex<double>, H5T_NATIVE_DOUBLE);
BOOSTSUB_H5_DATATYPE(std::complex<float>, H5T_NATIVE_FLOAT);

#else
typedef int hid_t;
typedef int herr_t;
typedef std::size_t hsize_t;
const int H5P_DEFAULT=0;

//return a non-sense integer
template <typename T>
inline hid_t
get_h5_datatype(const T&)
{
  return 0;
}

#endif
}
#endif

