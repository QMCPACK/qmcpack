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


#ifndef QMCPLUSPLUS_BOOSTADAPTOR_COLLECTIVES_OPERATIONS_H
#define QMCPLUSPLUS_BOOSTADAPTOR_COLLECTIVES_OPERATIONS_H

#include "type_traits/container_proxy.h"
#include "Message/mpi_datatype.h"

namespace qmcplusplus
{
namespace mpi
{
#if defined(HAVE_MPI)
//free function for allreduce
//   template<typename T, typename OP>
//   inline void all_reduce(const Communicate& comm,T& in, T& out)
//   {
//     container_proxy<T> t_in(in),t_out(out);
//     MPI_Datatype type_id=get_mpi_datatype(*t_in.data());
//     MPI_Allreduce(t_in.data(),t_out.data(),t_in.size(),type_id, OP, comm);
//   }

// /** generic function to perform allreduce */
// template<typename T, typename OP>
//   inline void all_reduce(const communicator& comm, T& in)
//   {
//     T out(in);
//     all_reduce<T,OP>(comm,in,out);
//     in=out;
//   }
//

template<typename T>
inline void reduce(const communicator& comm, T& in, T& out, int dest = 0)
{
  container_proxy<T> t_in(in), t_out(out);
  MPI_Datatype type_id = get_mpi_datatype(*t_in.data());
  MPI_Reduce(t_in.data(), t_out.data(), t_in.size(), type_id, MPI_SUM, dest, comm);
}

template<typename T>
inline void reduce(const communicator& comm, T& in, int dest = 0)
{
  T out(in);
  reduce<T>(comm, in, out, dest);
  in = out;
}


template<typename T, typename CT>
inline void all_gather(const communicator& comm, T& in, CT& out)
{
  container_proxy<T> t_in(in);
  container_proxy<CT> t_out(out);
  MPI_Datatype type_id = get_mpi_datatype(*t_in.data());
  int ierr             = MPI_Allgather(t_in.data(), t_in.size(), type_id, t_out.data(), t_in.size(), type_id, comm);
}

/** generic function to perform bcast
 *
 */


#else
template<typename T, typename OP>
inline void all_reduce(const communicator& comm, T& in, T& out)
{
  out = in;
}
template<typename T, typename OP>
inline void all_reduce(const communicator& comm, T& in)
{}
template<typename T, typename CT>
inline void all_gather(const communicator& comm, T& in, CT& out)
{
  out = in;
}



template<typename T>
inline void reduce(const communicator& comm, T& in, int dest = 0)
{}
template<typename T>
inline void reduce(const communicator& comm, T& in, T& out, int dest = 0)
{
  out = in;
}
#endif
} // namespace mpi
} // namespace qmcplusplus
#endif
