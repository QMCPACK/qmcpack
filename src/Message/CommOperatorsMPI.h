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
inline void Communicate::allreduce(T&)
{
  throw std::runtime_error("Need specialization for allreduce(T&)");
}

template<typename T>
inline void Communicate::reduce(T&)
{
  throw std::runtime_error("Need specialization for reduce(T&)");
}

template<typename T>
inline void Communicate::reduce(T* restrict, T* restrict, int n)
{
  throw std::runtime_error("Need specialization for reduce(T* restrict , T* restrict, int n)");
}

template<typename T>
inline void Communicate::reduce_in_place(T* restrict, int n)
{
  throw std::runtime_error("Need specialization for reduce_in_place(T* restrict, int n)");
}



template<typename T>
inline void Communicate::send(int dest, int tag, T&)
{
  throw std::runtime_error("Need specialization for send(int, int, T& )");
}


template<typename T>
inline void Communicate::allgather(T& sb, T& rb, int count)
{
  throw std::runtime_error("Need specialization for allgather(T&, T&, int)");
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
  MPI_Scatterv(t_in.data(), t_counts.data(), t_displ.data(), type_id, t_out.data(), t_out.size(), type_id, source, myMPI);
}

template<typename T>
inline Communicate::request Communicate::irecv(int source, int tag, T&)
{
  throw std::runtime_error("Need specialization for irecv(int source, int tag, T& )");
  return MPI_REQUEST_NULL;
}

template<typename T>
inline Communicate::request Communicate::isend(int dest, int tag, T&)
{
  throw std::runtime_error("Need specialization for isend(int source, int tag, T& )");
  return MPI_REQUEST_NULL;
}

template<typename T>
inline Communicate::request Communicate::irecv(int source, int tag, T*, int n)
{
  throw std::runtime_error("Need specialization for irecv(int source, int tag, T*, int )");
  return MPI_REQUEST_NULL;
}

template<typename T>
inline Communicate::request Communicate::isend(int dest, int tag, T*, int n)
{
  throw std::runtime_error("Need specialization for isend(int source, int tag, T*, int )");
  return MPI_REQUEST_NULL;
}

template<typename T>
inline void Communicate::allgather(T* sb, T* rb, int count)
{
  throw std::runtime_error("Need specialization for allgather(T*, T*, int)");
}

template<typename T, typename IT>
inline void Communicate::gatherv(T* sb, T* rb, int n, IT& counts, IT& displ, int dest)
{
  if (d_ncontexts == 1)
  {
    std::copy(sb, sb+n, rb);
    return;
  }
  qmcplusplus::container_proxy<IT> t_counts(counts), t_displ(displ);
  MPI_Datatype type_id = qmcplusplus::mpi::get_mpi_datatype(*sb);
  MPI_Gatherv(sb, n, type_id, rb, t_counts.data(), t_displ.data(), type_id, dest, myMPI);
}









template<>
inline void Communicate::allreduce(int& g)
{
  if (d_ncontexts == 1)
    return;
  int gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_INT, MPI_SUM, myMPI);
}

template<>
inline void Communicate::allreduce(long& g)
{
  if (d_ncontexts == 1)
    return;
  long gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_LONG, MPI_SUM, myMPI);
}

template<>
inline void Communicate::allreduce(unsigned long& g)
{
  if (d_ncontexts == 1)
    return;
  unsigned long gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_UNSIGNED_LONG, MPI_SUM, myMPI);
}

template<>
inline void Communicate::allreduce(float& g)
{
  if (d_ncontexts == 1)
    return;
  float gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_FLOAT, MPI_SUM, myMPI);
}

template<>
inline void Communicate::allreduce(double& g)
{
  if (d_ncontexts == 1)
    return;
  double gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_DOUBLE, MPI_SUM, myMPI);
}

template<>
inline void Communicate::allreduce(qmcplusplus::TinyVector<float, OHMMS_DIM>& g)
{
  if (d_ncontexts == 1)
    return;
  qmcplusplus::TinyVector<float, OHMMS_DIM> gt(g);
  MPI_Allreduce(g.begin(), gt.begin(), OHMMS_DIM, MPI_FLOAT, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(qmcplusplus::TinyVector<double, OHMMS_DIM>& g)
{
  if (d_ncontexts == 1)
    return;
  qmcplusplus::TinyVector<double, OHMMS_DIM> gt(g);
  MPI_Allreduce(g.begin(), gt.begin(), OHMMS_DIM, MPI_DOUBLE, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(qmcplusplus::TinyVector<int, OHMMS_DIM>& g)
{
  if (d_ncontexts == 1)
    return;
  qmcplusplus::TinyVector<int, OHMMS_DIM> gt(g);
  MPI_Allreduce(g.begin(), gt.begin(), OHMMS_DIM, MPI_INT, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(std::vector<int>& g)
{
  if (d_ncontexts == 1)
    return;
  std::vector<int> gt(g.size(), 0);
  MPI_Allreduce(g.data(), gt.data(), g.size(), MPI_INT, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(std::vector<long>& g)
{
  if (d_ncontexts == 1)
    return;
  std::vector<long> gt(g.size(), 0);
  MPI_Allreduce(g.data(), gt.data(), g.size(), MPI_LONG, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(std::vector<unsigned long>& g)
{
  if (d_ncontexts == 1)
    return;
  std::vector<unsigned long> gt(g.size(), 0);
  MPI_Allreduce(g.data(), gt.data(), g.size(), MPI_UNSIGNED_LONG, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(std::vector<float>& g)
{
  std::vector<float> gt(g.size(), 0.0f);
  MPI_Allreduce(g.data(), gt.data(), g.size(), MPI_FLOAT, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(std::vector<double>& g)
{
  std::vector<double> gt(g.size(), 0.0);
  MPI_Allreduce(g.data(), gt.data(), g.size(), MPI_DOUBLE, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(std::vector<std::complex<float>>& g)
{
  std::vector<std::complex<float>> gt(g.size(), std::complex<float>(0.0));
  MPI_Allreduce(g.data(), gt.data(), 2 * g.size(), MPI_FLOAT, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(std::vector<std::complex<double>>& g)
{
  std::vector<std::complex<double>> gt(g.size(), std::complex<double>(0.0));
  MPI_Allreduce(g.data(), gt.data(), 2 * g.size(), MPI_DOUBLE, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(PooledData<float>& g)
{
  PooledData<float> gt(g.size());
  MPI_Allreduce(g.data(), gt.data(), g.size(), MPI_FLOAT, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(PooledData<double>& g)
{
  PooledData<double> gt(g.size());
  MPI_Allreduce(g.data(), gt.data(), g.size(), MPI_DOUBLE, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void Communicate::allreduce(qmcplusplus::Matrix<float>& g)
{
  std::vector<float> gt(g.size());
  std::copy(g.begin(), g.end(), gt.begin());
  MPI_Allreduce(g.data(), gt.data(), g.size(), MPI_FLOAT, MPI_SUM, myMPI);
  std::copy(gt.begin(), gt.end(), g.data());
}

template<>
inline void Communicate::allreduce(qmcplusplus::Matrix<double>& g)
{
  std::vector<double> gt(g.size());
  copy(g.begin(), g.end(), gt.begin());
  MPI_Allreduce(g.data(), gt.data(), g.size(), MPI_DOUBLE, MPI_SUM, myMPI);
  copy(gt.begin(), gt.end(), g.data());
}

template<>
inline void Communicate::reduce(std::vector<float>& g)
{
  std::vector<float> gt(g.size(), 0.0f);
  MPI_Reduce(g.data(), gt.data(), g.size(), MPI_FLOAT, MPI_SUM, 0, myMPI);
  if (!d_mycontext)
    g = gt;
}

template<>
inline void Communicate::reduce(std::vector<double>& g)
{
  std::vector<double> gt(g.size(), 0.0);
  MPI_Reduce(g.data(), gt.data(), g.size(), MPI_DOUBLE, MPI_SUM, 0, myMPI);
  if (!d_mycontext)
    g = gt;
}

template<>
inline void Communicate::reduce(std::vector<int>& g)
{
  std::vector<int> gt(g.size(), 0.0);
  MPI_Reduce(g.data(), gt.data(), g.size(), MPI_INT, MPI_SUM, 0, myMPI);
  if (!d_mycontext)
    g = gt;
}

template<>
inline void Communicate::reduce(std::vector<long>& g)
{
  std::vector<long> gt(g.size(), 0.0);
  MPI_Reduce(g.data(), gt.data(), g.size(), MPI_LONG, MPI_SUM, 0, myMPI);
  if (!d_mycontext)
    g = gt;
}

template<>
inline void Communicate::reduce(int* restrict g, int* restrict res, int n)
{
  MPI_Reduce(g, res, n, MPI_INT, MPI_SUM, 0, myMPI);
}

template<>
inline void Communicate::reduce(double* restrict g, double* restrict res, int n)
{
  MPI_Reduce(g, res, n, MPI_DOUBLE, MPI_SUM, 0, myMPI);
}

template<>
inline void Communicate::reduce_in_place(double* restrict res, int n)
{
  if (!d_mycontext)
    MPI_Reduce(MPI_IN_PLACE, res, n, MPI_DOUBLE, MPI_SUM, 0, myMPI);
  else
    MPI_Reduce(res, NULL, n, MPI_DOUBLE, MPI_SUM, 0, myMPI);
}

template<>
inline void Communicate::reduce_in_place(float* restrict res, int n)
{
  if (!d_mycontext)
    MPI_Reduce(MPI_IN_PLACE, res, n, MPI_FLOAT, MPI_SUM, 0, myMPI);
  else
    MPI_Reduce(res, NULL, n, MPI_FLOAT, MPI_SUM, 0, myMPI);
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

template<>
inline void Communicate::send(int dest, int tag, std::vector<double>& g)
{
  MPI_Send(g.data(), g.size(), MPI_DOUBLE, dest, tag, myMPI);
}

template<>
inline Communicate::request Communicate::isend(int dest, int tag, std::vector<double>& g)
{
  request r;
  MPI_Isend(g.data(), g.size(), MPI_DOUBLE, dest, tag, myMPI, &r);
  return r;
}

template<>
inline Communicate::request Communicate::irecv(int source, int tag, std::vector<double>& g)
{
  request r;
  MPI_Irecv(g.data(), g.size(), MPI_DOUBLE, source, tag, myMPI, &r);
  return r;
}






template<>
inline void Communicate::allgather(std::vector<char>& sb, std::vector<char>& rb, int count)
{
  MPI_Allgather(sb.data(), count, MPI_CHAR, rb.data(), count, MPI_CHAR, myMPI);
}

template<>
inline void Communicate::allgather(std::vector<int>& sb, std::vector<int>& rb, int count)
{
  MPI_Allgather(sb.data(), count, MPI_INT, rb.data(), count, MPI_INT, myMPI);
}














template<>
inline void Communicate::allgather(char* sb, char* rb, int count)
{
  MPI_Allgather(sb, count, MPI_CHAR, rb, count, MPI_CHAR, myMPI);
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

template<>
inline void Communicate::allreduce(qmcplusplus::Matrix<std::complex<double>>& g)
{
  std::vector<std::complex<double>> gt(g.size());
  std::copy(g.begin(), g.end(), gt.begin());
  MPI_Allreduce(g.data(), gt.data(), 2 * g.size(), MPI_DOUBLE, MPI_SUM, myMPI);
  std::copy(gt.begin(), gt.end(), g.data());
}

template<>
inline void Communicate::allreduce(qmcplusplus::Matrix<std::complex<float>>& g)
{
  std::vector<std::complex<float>> gt(g.size());
  std::copy(g.begin(), g.end(), gt.begin());
  MPI_Allreduce(g.data(), gt.data(), 2 * g.size(), MPI_FLOAT, MPI_SUM, myMPI);
  std::copy(gt.begin(), gt.end(), g.data());
}



#endif
