//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_COMMUNICATION_OPERATORS_MPI_H
#define OHMMS_COMMUNICATION_OPERATORS_MPI_H
#include "Pools/PooledData.h"
#include "type_traits/container_proxy.h"
#include "Message/mpi_datatype.h"
#include <cstdint>
#include <stdexcept>
///dummy declarations to be specialized


template<typename T>
inline void Communicate::bcast(T& inout)
{
  if (d_ncontexts == 1)
    return;
  qmcplusplus::container_proxy<T> t_in(inout);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*t_in.data());
  MPI_Bcast(t_in.data(), t_in.size(), type_id, 0, myMPI);
}

template<typename T>
inline void Communicate::bcast(T* restrict inout, int n)
{
  if (d_ncontexts == 1)
    return;
  MPI_Bcast(inout, n, qmcplusplus::mpi::get_mpi_datatype(*inout), 0, myMPI);
}


template<typename T>
inline void Communicate::gather(T& sb, T& rb, int dest)
{
  if (d_ncontexts == 1)
  {
    rb = sb;
    return;
  }
  qmcplusplus::container_proxy<T> t_in(sb), t_out(rb);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*t_in.data());
  MPI_Gather(t_in.data(), t_in.size(), type_id, t_out.data(), t_in.size(), type_id, dest, myMPI);
}

template<typename T>
inline void Communicate::allreduce(T& g)
{
  if (d_ncontexts == 1)
    return;
  T gt(g);
  qmcplusplus::container_proxy<T> t_in(g), t_out(gt);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*t_in.data());
  MPI_Allreduce(t_in.data(), t_out.data(), t_in.size(), type_id, MPI_SUM, myMPI);
  g = gt;
}

template<typename T>
inline void Communicate::reduce(T& g)
{
  if (d_ncontexts == 1)
    return;
  T gt(g);
  qmcplusplus::container_proxy<T> t_in(g), t_out(gt);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*t_in.data());
  MPI_Reduce(t_in.data(), t_out.data(), t_in.size(), type_id, MPI_SUM, 0, myMPI);
  if (!d_mycontext)
    g = gt;
}


template<typename T>
inline void Communicate::reduce_in_place(T* restrict res, int n)
{
  if (d_ncontexts == 1)
    return;
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*res);
  if (!d_mycontext)
    MPI_Reduce(MPI_IN_PLACE, res, n, type_id, MPI_SUM, 0, myMPI);
  else
    MPI_Reduce(res, NULL, n, type_id, MPI_SUM, 0, myMPI);
}




template<typename T>
inline void Communicate::allgather(T& sb, T& rb, int count)
{
  if (d_ncontexts == 1)
  {
    rb = sb;
    return;
  }
  qmcplusplus::container_proxy<T> t_in(sb), t_out(rb);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*t_in.data());
  MPI_Allgather(t_in.data(), count, type_id, t_out.data(), count, type_id, myMPI);
}

template<typename T, typename IT>
inline void Communicate::gatherv(T& sb, T& rb, IT& counts, IT& displ, int dest)
{
  if (d_ncontexts == 1)
  {
    rb = sb;
    return;
  }
  qmcplusplus::container_proxy<T> t_in(sb), t_out(rb);
  qmcplusplus::container_proxy<IT> t_counts(counts), t_displ(displ);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*t_in.data());
  MPI_Gatherv(t_in.data(), t_in.size(), type_id, t_out.data(), t_counts.data(), t_displ.data(), type_id, dest, myMPI);
}

template<typename T>
inline void Communicate::scatter(T& sb, T& rb, int dest)
{
  if (d_ncontexts == 1)
  {
    rb = sb;
    return;
  }
  qmcplusplus::container_proxy<T> t_in(sb), t_out(rb);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*t_out.data());
  MPI_Scatter(t_in.data(), t_out.size(), type_id, t_out.data(), t_out.size(), type_id, dest, myMPI);
}

template<typename T, typename IT>
inline void Communicate::scatterv(T& sb, T& rb, IT& counts, IT& displ, int source)
{
  if (d_ncontexts == 1)
  {
    rb = sb;
    return;
  }
  qmcplusplus::container_proxy<T> t_in(sb), t_out(rb);
  qmcplusplus::container_proxy<IT> t_counts(counts), t_displ(displ);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*t_out.data());
  MPI_Scatterv(t_in.data(), t_counts.data(), t_displ.data(), type_id, t_out.data(), t_out.size(), type_id, source,
               myMPI);
}









template<typename T>
inline void Communicate::allgather(T* sb, T* rb, int count)
{
  if (d_ncontexts == 1)
  {
    for (int i = 0; i < count; ++i)
      rb[i] = sb[i];
    return;
  }
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*sb);
  MPI_Allgather(sb, count, type_id, rb, count, type_id, myMPI);
}

template<typename T, typename IT>
inline void Communicate::gatherv(T* sb, T* rb, int n, IT& counts, IT& displ, int dest)
{
  if (d_ncontexts == 1)
  {
    std::copy(sb, sb + n, rb);
    return;
  }
  qmcplusplus::container_proxy<IT> t_counts(counts), t_displ(displ);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*sb);
  MPI_Gatherv(sb, n, type_id, rb, t_counts.data(), t_displ.data(), type_id, dest, myMPI);
}


template<>
inline void Communicate::bcast(bool& g)
{
  int val = g ? 1 : 0;
  MPI_Bcast(&val, 1, MPI_INT, 0, myMPI);
  g = val != 0;
}


template<>
inline void Communicate::bcast(std::vector<bool>& g)
{
  std::vector<int> intVec(g.size());
  for (int i = 0; i < g.size(); i++)
    intVec[i] = g[i] ? 1 : 0;
  MPI_Bcast(&(intVec[0]), g.size(), MPI_INT, 0, myMPI);
  for (int i = 0; i < g.size(); i++)
    g[i] = intVec[i] != 0;
}


template<>
inline void Communicate::bcast(std::string& g)
{
  int string_size = g.size();

  bcast(string_size);
  if (rank() != 0)
    g.resize(string_size);

  bcast(g.data(), g.size());
}







template<typename T, typename TMPI, typename IT>
inline void Communicate::gatherv_in_place(T* buf, TMPI& datatype, IT& counts, IT& displ, int dest)
{
  if (!d_mycontext)
    MPI_Gatherv(MPI_IN_PLACE, 0, datatype, buf, counts.data(), displ.data(), datatype, dest, myMPI);
  else
    MPI_Gatherv(buf + displ[d_mycontext], counts[d_mycontext], datatype, NULL, counts.data(), displ.data(), datatype,
                dest, myMPI);
}


#endif
