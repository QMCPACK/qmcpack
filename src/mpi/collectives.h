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








/** generic function to perform bcast
 *
 */


#else








#endif
} // namespace mpi
} // namespace qmcplusplus
#endif
