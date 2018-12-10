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

#include <type_traits/container_proxy.h>
#include <mpi/mpi_datatype.h>
//#if defined(HAVE_MPI)
//#include <boost/mpi/operations.hpp>
//#endif

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
inline void reduce(const communicator& comm, T& in, T& out,  int dest=0)
{
  container_proxy<T> t_in(in),t_out(out);
  MPI_Datatype type_id=get_mpi_datatype(*t_in.data());
  MPI_Reduce(t_in.data(),t_out.data(),t_in.size(),type_id, MPI_SUM, dest, comm);
}

template<typename T>
inline void reduce(const communicator& comm, T& in, int dest=0)
{
  T out(in);
  reduce<T>(comm,in,out,dest);
  in=out;
}


/** generic function to perform allgather
 *
 * allgather of a scalar to a vectorized container
 */
template<typename T, typename CT>
inline void all_gather(const communicator& comm,T& in, CT& out)
{
  container_proxy<T> t_in(in);
  container_proxy<CT> t_out(out);
  MPI_Datatype type_id=get_mpi_datatype(*t_in.data());
  int ierr=MPI_Allgather(t_in.data(),t_in.size(),type_id,t_out.data(),t_in.size(),type_id,comm);
}

/** generic function to perform allgather
 *
 * allgather of a scalar to a vectorized container
 */
template<typename CT, typename IV>
inline void all_gatherv(const communicator& comm, CT& in, CT& out, IV& counts, IV& displ)
{
  container_proxy<CT> t_in(in),t_out(out);
  container_proxy<IV> t_counts(counts), t_displ(displ);
  MPI_Datatype type_id=get_mpi_datatype(*t_in.data());
  int ierr=MPI_Allgatherv(t_in.data(),t_in.size(),type_id,t_out.data(),t_counts.data(),t_displ.data(),type_id,comm);
}
/** generic function to perform allgather
 *
 * allgather of a scalar to a vectorized container
 */
template<typename CT, typename IV>
inline void gatherv(const communicator& comm, CT& in, CT& out, IV& counts, IV& displ, int dest=0)
{
  container_proxy<CT> t_in(in),t_out(out);
  container_proxy<IV> t_counts(counts), t_displ(displ);
  MPI_Datatype type_id=get_mpi_datatype(*t_in.data());
  int ierr=MPI_Gatherv(t_in.data(),t_in.size(),type_id,t_out.data(),t_counts.data(),t_displ.data(),type_id,dest,comm);
}

/** generic function to perform allgather
 *
 * allgather of a scalar to a vectorized container
 */
template<typename CT>
inline void gather(const communicator& comm, CT& in, CT& out, int dest=0)
{
  container_proxy<CT> t_in(in),t_out(out);
  MPI_Datatype type_id=get_mpi_datatype(*t_in.data());
  int ierr=MPI_Gather(t_in.data(),t_in.size(),type_id
                      ,t_out.data(),t_in.size(),type_id
                      ,dest,comm);
}

/** generic function to perform allgather
 *
 * allgather of a scalar to a vectorized container
 */
template<typename CT, typename IV>
inline void scatterv(const communicator& comm, CT& in, CT& out, IV& counts, IV& displ, int dest=0)
{
  container_proxy<CT> t_in(in),t_out(out);
  container_proxy<IV> t_counts(counts), t_displ(displ);
  MPI_Datatype type_id=get_mpi_datatype(*t_out.data());
  int ierr=MPI_Scatterv(t_in.data(),t_counts.data(),t_displ.data(),type_id
                        ,t_out.data(),t_out.size(),type_id,dest,comm);
}

/** generic function to perform allgather
 *
 * allgather of a scalar to a vectorized container
 */
template<typename CT>
inline void scatter(const communicator& comm, CT& in, CT& out, int dest=0)
{
  container_proxy<CT> t_in(in),t_out(out);
  MPI_Datatype type_id=get_mpi_datatype(*t_out.data());
  int ierr=MPI_Scatter(t_in.data(),t_out.size(),type_id
                       ,t_out.data(),t_out.size(),type_id,dest,comm);
}

/** generic function to perform bcast
 *
 */
template<typename CT>
inline void bcast(const communicator& comm, CT& inout, int source=0)
{
  if(comm.size()==1)
    return;
  container_proxy<CT> t_in(inout);
  MPI_Datatype type_id=get_mpi_datatype(*t_in.data());
  MPI_Bcast(t_in.data(),t_in.size(),type_id,source,comm);
}

template<typename T>
inline void bcast(const communicator& comm, T* inout, int n, int source=0)
{
  if(comm.size()==1)
    return;
  MPI_Bcast(inout,n,get_mpi_datatype(*inout),source,comm);
}

#else
template<typename T, typename OP>
inline void all_reduce(const communicator& comm,T& in, T& out)
{
  out=in;
}
template<typename T, typename OP>
inline void all_reduce(const communicator& comm,T& in) { }
template<typename T, typename CT>
inline void all_gather(const communicator& comm,T& in, CT& out)
{
  out=in;
}
template<typename CT, typename IT>
inline void all_gatherv(const communicator& comm, CT& in, CT& out, IT& counts, IT& displ)
{
  out=in;
}
template<typename CT, typename IT>
inline void gatherv(const communicator& comm, CT& in, CT& out, IT& counts, IT& displ, int dest=0)
{
  out=in;
}
template<typename CT>
inline void gather(const communicator& comm, CT& in, CT& out, int dest=0)
{
  out=in;
}
template<typename CT, typename IV>
inline void scatterv(const communicator& comm, CT& in, CT& out, IV& counts, IV& displ, int dest=0)
{
  out=in;
}
template<typename CT>
inline void scatter(const communicator& comm, CT& in, CT& out, int dest=0)
{
  out=in;
}
template<typename CT>
inline void bcast(const communicator& comm, CT& inout, int source=0)
{}
template<typename T>
inline void bcast(const communicator& comm, T* inout, int n, int source=0)
{}
template<typename T>
inline void reduce(const communicator& comm, T& in, int dest=0)
{}
template<typename T>
inline void reduce(const communicator& comm, T& in, T& out, int dest=0)
{
  out=in;
}
#endif
}
}
#endif
