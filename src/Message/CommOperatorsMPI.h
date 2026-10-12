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
#include "container_proxy.h"
#include <cstdint>
#include <stdexcept>
///dummy declarations to be specialized


template<typename T>
inline void Communicate::bcast(T& inout, int root)
{
  if (d_ncontexts == 1)
    return;
  qmcplusplus::container_proxy<T> t_in(inout);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*t_in.data());
  MPI_Bcast(t_in.data(), t_in.size(), type_id, root, myMPI);
}

template<typename T>
inline void Communicate::bcast(T* inout, int n, int root)
{
  if (d_ncontexts == 1)
    return;
  auto* addr = qmcplusplus::scalar_traits<T>::get_address(inout);
  MPI_Bcast(addr, n * qmcplusplus::scalar_traits<T>::DIM, qmcplusplus::mpi::get_mpi_datatype(*addr), root, myMPI);
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
inline void Communicate::reduce(T& g, int dest)
{
  if (d_ncontexts == 1)
    return;
  T gt(g);
  qmcplusplus::container_proxy<T> t_in(g), t_out(gt);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*t_in.data());
  MPI_Reduce(t_in.data(), t_out.data(), t_in.size(), type_id, MPI_SUM, dest, myMPI);
  if (d_mycontext == dest)
    g = gt;
}

template<typename T>
inline void Communicate::reduce(const T* sb, T* rb, int n, int dest)
{
  if (d_ncontexts == 1)
  {
    if (d_mycontext == dest)
      std::copy_n(sb, n, rb);
    return;
  }
  auto* s_addr         = qmcplusplus::scalar_traits<T>::get_address(sb);
  auto* r_addr         = qmcplusplus::scalar_traits<T>::get_address(rb);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*s_addr);
  MPI_Reduce(s_addr, r_addr, n * qmcplusplus::scalar_traits<T>::DIM, type_id, MPI_SUM, dest, myMPI);
}


template<typename T>
inline void Communicate::reduce_in_place(T* res, int n, int dest)
{
  if (d_ncontexts == 1)
    return;
  auto* addr           = qmcplusplus::scalar_traits<T>::get_address(res);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*addr);
  if (d_mycontext == dest)
    MPI_Reduce(MPI_IN_PLACE, addr, n * qmcplusplus::scalar_traits<T>::DIM, type_id, MPI_SUM, dest, myMPI);
  else
    MPI_Reduce(addr, NULL, n * qmcplusplus::scalar_traits<T>::DIM, type_id, MPI_SUM, dest, myMPI);
}


template<typename T>
inline void Communicate::allgather(T& sb, T& rb)
{
  if (d_ncontexts == 1)
  {
    rb = sb;
    return;
  }
  qmcplusplus::container_proxy<T> t_in(sb), t_out(rb);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*t_in.data());
  MPI_Allgather(t_in.data(), t_in.size(), type_id, t_out.data(), t_in.size(), type_id, myMPI);
}

template<typename T, typename IT>
inline void Communicate::gatherv(T& sb, T& rb, IT& counts, IT& displ, int dest)
{
  static_assert(qmcplusplus::scalar_traits<T>::DIM == 1, "Complex types not supported for this method");
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
  static_assert(qmcplusplus::scalar_traits<T>::DIM == 1, "Complex types not supported for this method");
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
  auto* addr_sb        = qmcplusplus::scalar_traits<T>::get_address(sb);
  auto* addr_rb        = qmcplusplus::scalar_traits<T>::get_address(rb);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*addr_sb);
  MPI_Allgather(addr_sb, count * qmcplusplus::scalar_traits<T>::DIM, type_id, addr_rb,
                count * qmcplusplus::scalar_traits<T>::DIM, type_id, myMPI);
}

template<typename T, typename IT>
inline void Communicate::gatherv(T* sb, T* rb, int n, IT& counts, IT& displ, int dest)
{
  static_assert(qmcplusplus::scalar_traits<T>::DIM == 1, "Complex types not supported for this method");
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
inline void Communicate::bcast(bool& g, int root)
{
  int val = g ? 1 : 0;
  MPI_Bcast(&val, 1, MPI_INT, root, myMPI);
  g = val != 0;
}


template<>
inline void Communicate::bcast(std::vector<bool>& g, int root)
{
  std::vector<int> intVec(g.size());
  for (int i = 0; i < g.size(); i++)
    intVec[i] = g[i] ? 1 : 0;
  MPI_Bcast(&(intVec[0]), g.size(), MPI_INT, root, myMPI);
  for (int i = 0; i < g.size(); i++)
    g[i] = intVec[i] != 0;
}


template<>
inline void Communicate::bcast(std::string& g, int root)
{
  int string_size = g.size();

  bcast(string_size, root);
  if (rank() != root)
    g.resize(string_size);

  bcast(g.data(), g.size(), root);
}


template<typename T, typename TMPI, typename IT>
inline void Communicate::gatherv_in_place(T* buf, const TMPI& datatype, IT& counts, IT& displ, int dest)
{
  static_assert(qmcplusplus::scalar_traits<T>::DIM == 1, "Complex types not supported for this method");
  if (!d_mycontext)
    MPI_Gatherv(MPI_IN_PLACE, 0, datatype, buf, counts.data(), displ.data(), datatype, dest, myMPI);
  else
    MPI_Gatherv(buf + displ[d_mycontext], counts[d_mycontext], datatype, NULL, counts.data(), displ.data(), datatype,
                dest, myMPI);
}


#endif
