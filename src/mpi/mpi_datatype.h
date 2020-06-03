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


#ifndef QMCPLUSPLUS_MPI_DATATYPEDEFINE_H
#define QMCPLUSPLUS_MPI_DATATYPEDEFINE_H


#include <type_traits/scalar_traits.h>
#if defined(HAVE_MPI)
#include <mpi.h>
#else
typedef int MPI_Datatype;
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
typedef MPI_Status status;

template<typename T>
inline MPI_Datatype get_mpi_datatype(const T&)
{
  return MPI_BYTE;
}

#define BOOSTSUB_MPI_DATATYPE(CppType, MPITYPE)                 \
  template<>                                                    \
  inline MPI_Datatype get_mpi_datatype<CppType>(const CppType&) \
  {                                                             \
    return MPITYPE;                                             \
  }

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

template<typename T>
void free_column_type(T& datatype)
{
  MPI_Type_free(&datatype);
}

template<typename T>
MPI_Datatype construct_column_type(const T* element, int nrow, int ncol)
{
  MPI_Datatype column_type;
  MPI_Datatype column_type1;
  MPI_Type_vector(nrow, 1, ncol, get_mpi_datatype(*element), &column_type);
  MPI_Type_commit(&column_type);
  MPI_Type_create_resized(column_type, 0, sizeof(T), &column_type1);
  MPI_Type_commit(&column_type1);
  free_column_type(column_type);
  return column_type1;
}

#else
typedef int status;
typedef int request;
typedef int MPI_Datatype;

//return a non-sense integer
template<typename T>
inline MPI_Datatype get_mpi_datatype(const T&)
{
  return 0;
}

template<typename T>
void free_column_type(T& datatype)
{}

template<typename T>
MPI_Datatype construct_column_type(const T* element, int nrow, int ncol)
{
  return 0;
}

#endif
} // namespace mpi
} // namespace qmcplusplus
#endif
