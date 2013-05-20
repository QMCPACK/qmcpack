//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file  mpi_datatype.h
 *
 * data type compatible with MPI. Provides minimalistic interfaces to MPI
 * Do not use "using namespace mpi". Always use with mpi::
 * Expect to be replaced by boost::mpi
 */
#ifndef QMCPLUSPLUS_MPI_DATATYPEDEFINE_H
#define QMCPLUSPLUS_MPI_DATATYPEDEFINE_H


#include <type_traits/scalar_traits.h>
#if defined(HAVE_MPI)
#include <mpi.h>
#endif

namespace qmcplusplus
{
namespace mpi
{

typedef Communicate communicator;

#if defined(HAVE_MPI)

///@typedef mpi::request
typedef MPI_Request request;
///@typedef mpi::status
typedef MPI_Status  status;

template <typename T>
inline MPI_Datatype
get_mpi_datatype(const T&)
{
  return MPI_BYTE;
}

#define BOOSTSUB_MPI_DATATYPE(CppType, MPITYPE)                \
template<>                                                     \
inline MPI_Datatype                                             \
get_mpi_datatype< CppType >(const CppType&) { return MPITYPE; }

BOOSTSUB_MPI_DATATYPE(short, MPI_SHORT);

BOOSTSUB_MPI_DATATYPE(int, MPI_INT);

BOOSTSUB_MPI_DATATYPE(long, MPI_LONG);

BOOSTSUB_MPI_DATATYPE(float, MPI_FLOAT);

BOOSTSUB_MPI_DATATYPE(double, MPI_DOUBLE);

BOOSTSUB_MPI_DATATYPE(long double, MPI_LONG_DOUBLE);

BOOSTSUB_MPI_DATATYPE(unsigned char, MPI_UNSIGNED_CHAR);

BOOSTSUB_MPI_DATATYPE(unsigned short, MPI_UNSIGNED_SHORT);

BOOSTSUB_MPI_DATATYPE(unsigned int, MPI_UNSIGNED);

BOOSTSUB_MPI_DATATYPE(unsigned long, MPI_UNSIGNED_LONG);

BOOSTSUB_MPI_DATATYPE(std::complex<double>, MPI_DOUBLE);

BOOSTSUB_MPI_DATATYPE(std::complex<float>, MPI_FLOAT);

#else
typedef int  status;
typedef int  request;
typedef int MPI_Datatype;
//return a non-sense integer
template <typename T>
inline MPI_Datatype get_mpi_datatype(const T&)
{
  return 0;
}
#endif
}
}//end of mpi
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 894 $   $Date: 2006-02-03 10:52:38 -0600 (Fri, 03 Feb 2006) $
 * $Id: hdf_datatype.h 894 2006-02-03 16:52:38Z jnkim $
 ***************************************************************************/

