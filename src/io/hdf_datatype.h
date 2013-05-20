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
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 894 $   $Date: 2006-02-03 10:52:38 -0600 (Fri, 03 Feb 2006) $
 * $Id: hdf_datatype.h 894 2006-02-03 16:52:38Z jnkim $
 ***************************************************************************/

