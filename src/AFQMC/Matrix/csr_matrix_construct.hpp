////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef SPARSE_CSR_MATRIX_CONSTRUCT_HPP
#define SPARSE_CSR_MATRIX_CONSTRUCT_HPP

#include <complex>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
#include <tuple>
#include <numeric>
#include <memory>
#include <type_traits> // enable_if
#include <algorithm>
#include <mutex>

#include "Configuration.h"
#include "Utilities/FairDivide.h"
#include "AFQMC/config.h"
#include "AFQMC/Utilities/Utils.hpp"
#include "mpi3/shared_communicator.hpp"
#include "mpi3/shared_window.hpp"
#include "mpi3/shm/mutex.hpp"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Utilities/myTimer.h"
#include "AFQMC/Numerics/detail/utilities.hpp"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"

namespace csr
{
namespace shm
{
/*
 * Constructs a csr_matrix from the elements of the multi::array M applying a truncation.
 * The input matrix is only meaningful in the core with rank 0
 */
template<class CSR,
         class MultiArray2D,
         typename = typename std::enable_if<std::decay<MultiArray2D>::type::dimensionality == 2>::type>
CSR construct_csr_matrix_single_input(MultiArray2D&& M, double cutoff, char TA, boost::mpi3::shared_communicator& comm)
{
  assert(TA == 'N' || TA == 'T' || TA == 'H');
  std::vector<std::size_t> counts;
  using int_type = typename CSR::index_type;
  int_type nr, nc;
  if (comm.rank() == 0)
  {
    if (TA == 'N')
    {
      nr = std::get<0>(M.sizes());
      nc = std::get<1>(M.sizes());
      counts.resize(nr);
      for (int_type i = 0; i < nr; i++)
        for (int_type j = 0; j < nc; j++)
          if (std::abs(M[i][j]) > cutoff)
            ++counts[i];
    }
    else
    {
      nr = std::get<1>(M.sizes());
      nc = std::get<0>(M.sizes());
      counts.resize(nr);
      for (int_type i = 0; i < std::get<0>(M.sizes()); i++)
        for (int_type j = 0; j < std::get<1>(M.sizes()); j++)
          if (std::abs(M[i][j]) > cutoff)
            ++counts[j];
    }
  }
  comm.broadcast_value(nr);
  comm.broadcast_value(nc);
  if (comm.rank() != 0)
    counts.resize(nr);
  comm.broadcast_n(counts.begin(), counts.size());

  CSR csr_mat(std::tuple<std::size_t, std::size_t>{nr, nc}, std::tuple<std::size_t, std::size_t>{0, 0}, counts,
              qmcplusplus::afqmc::shared_allocator<typename CSR::value_type>(comm));

  if (comm.rank() == 0)
  {
    if (TA == 'N')
    {
      for (int_type i = 0; i < nr; i++)
        for (int_type j = 0; j < nc; j++)
          if (std::abs(M[i][j]) > cutoff)
            csr_mat.emplace_back({i, j}, static_cast<typename CSR::value_type>(M[i][j]));
    }
    else if (TA == 'T')
    {
      for (int_type i = 0; i < std::get<1>(M.sizes()); i++)
        for (int_type j = 0; j < std::get<0>(M.sizes()); j++)
          if (std::abs(M[j][i]) > cutoff)
            csr_mat.emplace_back({i, j}, static_cast<typename CSR::value_type>(M[j][i]));
    }
    else if (TA == 'H')
    {
      for (int_type i = 0; i < std::get<1>(M.sizes()); i++)
        for (int_type j = 0; j < std::get<0>(M.sizes()); j++)
          if (std::abs(M[j][i]) > cutoff)
            csr_mat.emplace_back({i, j}, static_cast<typename CSR::value_type>(ma::conj(M[j][i])));
    }
  }
  csr_mat.remove_empty_spaces();
  comm.barrier();

  return csr_mat;
}

/*
 * Constructs a csr_matrix from the elements in the container (of tuples).  
 * All cores in the communicator contribute values. 
 */
template<class CSR, class Container>
CSR construct_csr_matrix_multiple_input(Container const& M,
                                        std::size_t nr,
                                        std::size_t nc,
                                        char TA,
                                        boost::mpi3::shared_communicator& comm)
{
  myTimer Timer;
  Timer.reset("G0");
  Timer.start("G0");

  assert(TA == 'N' || TA == 'T' || TA == 'H');
  using std::get;
  using VType = typename CSR::value_type;
  using IType = typename CSR::index_type;
  using PType = typename CSR::int_type;
  using UCSR  = typename CSR::base;

  std::vector<std::size_t> counts(nr);
  if (TA == 'N')
    for (auto& v : M)
      ++counts[get<0>(v)];
  else
    for (auto& v : M)
      ++counts[get<1>(v)];
  comm.all_reduce_in_place_n(counts.begin(), counts.size(), std::plus<>());

  Timer.stop("G0");
  qmcplusplus::app_log() << " In construct_csr_matrix_multiple_input allreduce: " << Timer.total("G0") << std::endl;
  Timer.reset("G0");
  Timer.start("G0");

  UCSR ucsr_mat(std::tuple<std::size_t, std::size_t>{nr, nc}, std::tuple<std::size_t, std::size_t>{0, 0}, counts,
                qmcplusplus::afqmc::shared_allocator<VType>(comm));

  for (std::size_t r = 0; r < comm.size(); r++)
  {
    comm.barrier();
    if (comm.rank() == r)
    {
      if (TA == 'N')
        for (auto& v : M)
          ucsr_mat.emplace({get<0>(v), get<1>(v)}, get<2>(v));
      else if (TA == 'T')
        for (auto& v : M)
          ucsr_mat.emplace({get<1>(v), get<0>(v)}, get<2>(v));
      else if (TA == 'H')
        for (auto& v : M)
          ucsr_mat.emplace({get<1>(v), get<0>(v)}, ma::conj(get<2>(v)));
    }
    comm.barrier();
  }

  Timer.stop("G0");
  qmcplusplus::app_log() << " In construct_csr_matrix_multiple_input emplace: " << Timer.total("G0") << std::endl;

  return CSR(std::move(ucsr_mat));
}

/*
 * Constructs a new csr_matrix by adding to the given csr_matrix::base all the elements
 * on the other nodes in TG.Global().  
 * All nodes (including all TGs) contribute values. 
 * All nodes involved in the calculation obtain identical copies of the resulting csr matrix. 
 */
template<class CSR, class task_group>
CSR construct_csr_matrix_from_distributed_ucsr(typename CSR::base&& ucsr, task_group& TG)
{
  if (ucsr.size(0) == 0 || ucsr.size(1) == 0)
    return CSR(std::move(ucsr));
  ;
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();
  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();

  std::size_t ak0, ak1;
  std::tie(ak0, ak1) = FairDivideBoundary(std::size_t(coreid), ucsr.size(0), std::size_t(ncores));

  std::vector<std::size_t> counts_local(ak1 - ak0);
  std::vector<std::size_t> counts_global(ucsr.size(0));
  for (std::size_t r = ak0, r0 = 0; r < ak1; ++r, ++r0)
  {
    counts_local[r0] = std::size_t(*ucsr.pointers_end(r) - *ucsr.pointers_begin(r));
    counts_global[r] = counts_local[r0];
  }
  TG.Global().all_reduce_in_place_n(counts_global.begin(), counts_global.size(), std::plus<>());

  if (std::accumulate(counts_global.begin(), counts_global.end(), std::size_t(0)) == 0)
    return CSR(std::move(ucsr));
  ;

  // if this has been previously "reserved", it will do nothing
  ucsr.reserve(counts_global);

  std::size_t nterms = std::accumulate(counts_local.begin(), counts_local.end(), std::size_t(0));
  auto sz_per_node   = TG.Cores().all_gather_value(nterms);
  std::size_t maxnterms(*std::max_element(sz_per_node.begin(), sz_per_node.end()));
  std::vector<typename CSR::value_type> vals(maxnterms);
  std::vector<typename CSR::index_type> cols(maxnterms);
  std::vector<std::size_t> counts(counts_local.size());
  for (int ni = 0; ni < nnodes; ++ni)
  {
    if (nodeid == ni)
    {
      auto v0 = vals.begin();
      auto c0 = cols.begin();
      std::copy(counts_local.begin(), counts_local.end(), counts.begin());
      for (std::size_t r = ak0, r0 = 0; r < ak1; ++r, ++r0)
      {
        v0 = std::copy_n(to_address(ucsr.non_zero_values_data(r)), counts[r0], v0);
        c0 = std::copy_n(to_address(ucsr.non_zero_indices2_data(r)), counts[r0], c0);
      }
    }
    TG.Cores().broadcast_n(vals.begin(), sz_per_node[ni], ni);
    TG.Cores().broadcast_n(cols.begin(), sz_per_node[ni], ni);
    TG.Cores().broadcast_n(counts.begin(), counts.size(), ni);
    if (nodeid != ni)
    {
      auto itC                   = cols.begin();
      auto itV                   = vals.begin();
      typename CSR::index_type r = static_cast<typename CSR::index_type>(ak0);
      for (auto i : counts)
      {
        while (i-- > 0)
          ucsr.emplace({r, *itC++}, *itV++);
        ++r;
      }
    }
  }
  TG.node_barrier();
  return CSR(std::move(ucsr));
}

/*
 * Constructs a new csr_matrix from the elements in the container Q. 
 * All cores (including all TGs) contribute values. 
 * All nodes involved in the calculation obtain identical copies of the resulting csr matrix. 
 */
template<class CSR, class Container, class task_group>
CSR construct_csr_matrix_from_distributed_containers(Container const& Q,
                                                     std::size_t nr,
                                                     std::size_t nc,
                                                     task_group& TG,
                                                     bool needs_lock = true)
{
  using std::get;
  using Type = typename Container::value_type;
  if (nr == 0 || nc == 0)
    return CSR(std::tuple<std::size_t, std::size_t>{nr, nc}, std::tuple<std::size_t, std::size_t>{0, 0}, 0,
               typename CSR::alloc_type(TG.Node()));
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();
  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();

  std::vector<std::size_t> counts_global(nr);
  for (auto& v : Q)
    ++counts_global[get<0>(v)];
  TG.Global().all_reduce_in_place_n(counts_global.begin(), counts_global.size(), std::plus<>());

  if (std::accumulate(counts_global.begin(), counts_global.end(), std::size_t(0)) == 0)
    return CSR(std::tuple<std::size_t, std::size_t>{nr, nc}, std::tuple<std::size_t, std::size_t>{0, 0}, 0,
               typename CSR::alloc_type(TG.Node()));

  typename CSR::base ucsr(std::tuple<std::size_t, std::size_t>{nr, nc}, std::tuple<std::size_t, std::size_t>{0, 0},
                          counts_global, typename CSR::alloc_type(TG.Node()));

  std::size_t nterms = Q.size();
  auto sz_per_node   = TG.Cores().all_gather_value(nterms);
  std::size_t maxnterms(*std::max_element(sz_per_node.begin(), sz_per_node.end()));
  std::vector<Type> Qc;
  Qc.reserve(maxnterms);
  boost::mpi3::shm::mutex m(TG.Node());
  for (int ni = 0; ni < nnodes; ++ni)
  {
    Qc.resize(sz_per_node[ni]);
    if (nodeid == ni)
      std::copy_n(Q.begin(), sz_per_node[ni], Qc.begin());
    //TG.Cores().broadcast_n(Qc.begin(),sz_per_node[ni],ni);
    // replace when tuples can me communicated directly
    MPI_Bcast(Qc.data(), sz_per_node[ni] * sizeof(Type), MPI_CHAR, ni, &(TG.Cores()));
    if (needs_lock)
    {
      std::lock_guard<boost::mpi3::shm::mutex> guard(m);
      for (auto& v : Qc)
        ucsr.emplace({get<0>(v), get<1>(v)}, get<2>(v));
    }
    else
    {
      for (auto& v : Qc)
        ucsr.emplace({get<0>(v), get<1>(v)}, get<2>(v));
    }
  }
  TG.node_barrier();
  return CSR(std::move(ucsr));
}

/*
 * Constructs a new csr_matrix from the elements in the container Q. 
 * The global matrix (including all elements in all cores) will be evenly distributed
 * across the nodes in every task group. No particular structure will be followed in the
 * partitioning, only strict distribution of non-zero elements.
 * All TGs will have identical distributions among its nodes. 
 * This approach uses more memory (up to 2 copies of the submatrix), but avoids
 * having to execute code to count the number of terms ahead of time.
 * Current algorithm:
 *   1. all_gather on equivalent nodes on a TG. This will leave every TG with the entire matrix.
 *   2. Count elements and determine final nnz per node.
 *   3. Exchange elements between nodes to achieve desired distribution
 *        i. Cores with too many elements sends excess to head cores of nodes with too few   
 *   4. construct CSR   
 */
template<class CSR, class Container, class task_group>
CSR construct_distributed_csr_matrix_from_distributed_containers(Container& Q,
                                                                 std::size_t nr,
                                                                 std::size_t nc,
                                                                 task_group& TG)
{
  myTimer Timer;
  Timer.reset("G0");
  Timer.start("G0");

  using std::get;
  using Type = typename Container::value_type;
  if (nr == 0 || nc == 0)
    return CSR(std::tuple<std::size_t, std::size_t>{nr, nc}, std::tuple<std::size_t, std::size_t>{0, 0}, 0,
               typename CSR::alloc_type(TG.Node()));
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();
  int nodeid        = TG.getNodeID();
  int node_number   = TG.getLocalGroupNumber();
  int nnodes_per_TG = TG.getNGroupsPerTG();

  // 1. Define new communicator for equivalent cores
  boost::mpi3::communicator eq_cores(TG.Cores().split(node_number, TG.Cores().rank()));
  // this comm is similar to TG.TG, but includes all the cores in the node in the same comm
  // it is used to balance the distribution within the nodes on each TG.TG
  boost::mpi3::communicator eq_node_group(TG.Global().split(nodeid / nnodes_per_TG, TG.Global().rank()));

  // 2. count elements
  long nterms_total            = (TG.Global() += Q.size());
  auto nterms_per_core         = eq_cores.all_gather_value(Q.size());
  long nterms_in_local_core    = std::accumulate(nterms_per_core.begin(), nterms_per_core.end(), long(0));
  auto nterms_per_core_in_node = TG.Node().all_gather_value(nterms_in_local_core);
  long nterms_node_before_loadbalance =
      std::accumulate(nterms_per_core_in_node.begin(), nterms_per_core_in_node.end(), long(0));
  std::vector<long> bounds(nnodes_per_TG + 1);
  FairDivide(size_t(nterms_total), nnodes_per_TG, bounds);

  // all extra terms are received by CoreID==0
  if (coreid == 0)
    Q.reserve(nterms_in_local_core +
              std::max(long(0), bounds[node_number + 1] - bounds[node_number] - nterms_node_before_loadbalance));
  else
    Q.reserve(nterms_in_local_core);

  // 3. all_gather_v on eq_cores (use mpi3 when tested and ready)
  //     i. move data in Q to allow use of MPI_IN_PLACE in Allgatherv
  Q.resize(nterms_in_local_core);
  if (eq_cores.rank() > 0)
  {
    size_t n = std::accumulate(nterms_per_core.begin(), nterms_per_core.begin() + eq_cores.rank() + 1, size_t(0));
    std::copy_backward(Q.begin(), Q.begin() + nterms_per_core[eq_cores.rank()], Q.begin() + n);
  }

  Timer.stop("G0");
  qmcplusplus::app_log() << " In construct_distributed_csr_matrix_from_distributed_containers setup: "
                         << Timer.total("G0") << std::endl;
  Timer.reset("G0");
  Timer.start("G0");

  std::vector<int> recvcounts(eq_cores.size());
  std::vector<int> displ(eq_cores.size());
  long cnt(0);
  for (int p = 0; p < eq_cores.size(); p++)
  {
    if (nterms_per_core[p] * sizeof(Type) >= static_cast<long>(std::numeric_limits<int>::max()))
      throw std::out_of_range("row size exceeded the maximum");
    recvcounts[p] = int(nterms_per_core[p] * sizeof(Type));
    displ[p]      = int(cnt * sizeof(Type));
    cnt += nterms_per_core[p];
    if (cnt * sizeof(Type) >= static_cast<long>(std::numeric_limits<int>::max()))
      throw std::out_of_range("row size exceeded the maximum");
  }
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, Q.data(), recvcounts.data(), displ.data(), MPI_CHAR, eq_cores.get());

  Timer.stop("G0");
  qmcplusplus::app_log() << " In construct_distributed_csr_matrix_from_distributed_containers Allgatherv: "
                         << Timer.total("G0") << std::endl;
  Timer.reset("G0");
  Timer.start("G0");

  boost::multi::array<long, 2> counts_per_core({nnodes_per_TG, ncores});
  std::fill_n(counts_per_core.origin(), ncores * nnodes_per_TG, long(0));
  counts_per_core[node_number][coreid] = nterms_in_local_core;
  eq_node_group.all_reduce_in_place_n(counts_per_core.origin(), ncores * nnodes_per_TG, std::plus<>());

  // calculate how many terms each core is potentially sending
  // when algorithm works, make nplus a vector of pairs to save space, no need to store mostly zeros
  boost::multi::array<long, 2> nplus({nnodes_per_TG, ncores});
  std::fill_n(nplus.origin(), ncores * nnodes_per_TG, long(0));
  auto tgrank = eq_node_group.rank();
  // number of terms I'm sending
  long deltaN = 0;
  // excess terms before me
  long nbefore = 0;
  for (int ni = 0, ip = 0; ni < nnodes_per_TG; ++ni)
  {
    long nuptome = 0;
    for (int ci = 0; ci < ncores; ++ci, ++ip)
    {
      nuptome += counts_per_core[ni][ci];
      long dn = std::min(counts_per_core[ni][ci], nuptome - (bounds[ni + 1] - bounds[ni]));
      if (dn > long(0))
      {
        if (ip < tgrank)
          nbefore += long(dn);
        if (ip == tgrank)
          deltaN = long(dn);
        nplus[ni][ci] = dn;
      }
    }
  }

  // calculate how many terms each head core is expecting
  std::vector<long> nminus(nnodes_per_TG);
  for (int i = 0; i < nnodes_per_TG; i++)
  {
    long dn = (bounds[i + 1] - bounds[i]) -
        std::accumulate(counts_per_core[i].origin(), counts_per_core[i].origin() + ncores, long(0));
    if (dn > 0)
      nminus[i] = long(dn);
  }

  Timer.stop("G0");
  qmcplusplus::app_log() << " In construct_distributed_csr_matrix_from_distributed_containers setup comm: "
                         << Timer.total("G0") << std::endl;
  Timer.reset("G0");
  Timer.start("G0");

  /*
qmcplusplus::app_log()<<nnodes_per_TG <<" " <<ncores <<"\n";
for(int i=0; i<nnodes_per_TG; i++)
 qmcplusplus::app_log()<<nminus[i] <<" ";
qmcplusplus::app_log()<<std::endl;
for(int i=0; i<nnodes_per_TG; i++) {
 for(int j=0; j<ncores; j++)
  qmcplusplus::app_log()<<nplus[i][j] <<" ";
 qmcplusplus::app_log()<<std::endl;
}

std::ofstream out((std::string("debug.")+std::to_string(TG.Global().rank())).c_str());
*/
  //
  if (deltaN > long(0))
  {
    // I'm sending somewhere
    long nsent(0);
    for (int i = 0; i < nnodes_per_TG; i++)
    {
      if (nminus[i] == long(0))
        continue;
      if (nsent + nminus[i] > nbefore)
      {
        // send to i
        long to_send = std::min(deltaN, nsent + nminus[i] - nbefore);
        if (to_send * sizeof(Type) >= static_cast<long>(std::numeric_limits<int>::max()))
          throw std::out_of_range("row size exceeded the maximum");
        int to_send_int = int(to_send * sizeof(Type));
        //out<<" sending " <<to_send <<" " <<deltaN <<" " <<i*ncores <<" " <<i <<std::endl;
        MPI_Send(Q.data() + Q.size() - to_send, to_send_int, MPI_CHAR, i * ncores, tgrank, eq_node_group.get());
        //out<<" done sending " <<to_send <<" " <<deltaN <<" " <<i*ncores <<" " <<i <<std::endl;
        nbefore += to_send;
        deltaN -= to_send;
        while (to_send-- > long(0))
          Q.pop_back();
      }
      nsent += nminus[i];
      if (deltaN == long(0))
        break;
    }
    if (deltaN > long(0))
      throw std::out_of_range("detlaN > 0");
  }
  else if (coreid == 0 && nminus[node_number] > long(0))
  {
    deltaN  = nminus[node_number];
    nbefore = std::accumulate(nminus.begin(), nminus.begin() + node_number, long(0));
    long nsent(0);
    MPI_Status st;
    for (int ni = 0, ip = 0; ni < nnodes_per_TG; ni++)
    {
      for (int ci = 0; ci < ncores; ++ci, ++ip)
      {
        if (nplus[ni][ci] == long(0))
          continue;
        if (nsent + nplus[ni][ci] > nbefore)
        {
          // receive from ip
          long to_recv = std::min(deltaN, nsent + nplus[ni][ci] - nbefore);
          if (to_recv * sizeof(Type) >= static_cast<long>(std::numeric_limits<int>::max()))
            throw std::out_of_range("row size exceeded the maximum");
          int to_recv_int = int(to_recv * sizeof(Type));
          long curr_sz    = Q.size();
          Q.resize(curr_sz + to_recv);
          //out<<" receiving " <<to_recv <<" " <<deltaN <<" " <<nbefore <<" " <<ip <<std::endl;
          MPI_Recv(Q.data() + curr_sz, to_recv_int, MPI_CHAR, ip, ip, eq_node_group.get(), &st);
          nbefore += to_recv;
          deltaN -= to_recv;
          //out<<" done receiving " <<to_recv <<" " <<deltaN <<" " <<nbefore <<" " <<ip <<std::endl;
        }
        nsent += nplus[ni][ci];
        if (deltaN == long(0))
          break;
      }
      if (deltaN == long(0))
        break;
    }
    if (deltaN > long(0))
      throw std::out_of_range("detlaN > 0");
  }
  else
  {
    //out<<" nothing to do. \n";
  }

  //out.close();

  Timer.stop("G0");
  qmcplusplus::app_log() << " In construct_distributed_csr_matrix_from_distributed_containers send/recv: "
                         << Timer.total("G0") << std::endl;
  Timer.reset("G0");
  Timer.start("G0");


  long final_nterms_node = (TG.Node() += Q.size());
  std::vector<long> final_counts(TG.Cores().size());
  TG.Cores().gather_n(&final_nterms_node, 1, final_counts.data(), 0);
  if (TG.Global().root())
  {
    qmcplusplus::app_log() << " In construct_distributed_csr_matrix_from_distributed_containers: \n";
    qmcplusplus::app_log() << " Partitioning of CSR matrix over TG (nnz per node): ";
    for (int i = 0; i < nnodes_per_TG; i++)
      qmcplusplus::app_log() << final_counts[i] << " ";
    qmcplusplus::app_log() << std::endl;
  }

  TG.node_barrier();
  return construct_csr_matrix_multiple_input<CSR>(Q, nr, nc, 'N', TG.Node());
}

// assume partition by row
template<class CSR, class task_group, typename integer = int>
typename CSR::template matrix_view<integer> local_balanced_partition(CSR& M, task_group& TG)
{
  using std::size_t;
  using array_ = std::array<size_t, 4>;
  if (TG.getNCoresPerTG() == 1)
  {
    return M[array_{0, M.size(0), 0, M.size(1)}];
  }
  else
  {
    std::vector<size_t> bounds(TG.getNCoresPerTG() + 1);
    if (TG.getCoreID() == 0)
    {
      std::vector<size_t> indx(M.size(0) + 1);
      indx[0] = size_t(0);
      for (size_t i = 0, cnt = 0; i < M.size(0); i++)
      {
        cnt += size_t(M.pointers_end()[i] - M.pointers_begin()[i]);
        indx[i + 1] = cnt;
      }
      qmcplusplus::afqmc::balance_partition_ordered_set(indx, bounds);
      if (TG.Global().root())
      {
        qmcplusplus::app_log() << std::endl << "Partition of csr matrix (rows) over cores in local TG: ";
        for (int i = 0; i <= TG.getNCoresPerTG(); i++)
          qmcplusplus::app_log() << bounds[i] << " ";
        qmcplusplus::app_log() << std::endl;
        qmcplusplus::app_log() << "Number of terms in csr matrix view per core: ";
        for (int i = 0; i < TG.getNCoresPerTG(); i++)
          qmcplusplus::app_log() << indx[bounds[i + 1]] - indx[bounds[i]] << " ";
        qmcplusplus::app_log() << std::endl << std::endl;
      }
    }
    TG.Node().broadcast_n(bounds.begin(), bounds.size());
    return M[array_{bounds[TG.getLocalTGRank()], bounds[TG.getLocalTGRank() + 1], 0, M.size(1)}];
    //return M[array_{0,M.size(0),0,M.size(1)}];
  }
}

/*
 * Constructs a vector of csr_matrix as a copy from a given csr_matrix but casted to single precision
 */
template<class CSRsp, class CSR>
std::vector<CSRsp> CSRvector_to_single_precision(std::vector<CSR> const& A)
{
  using Alloc = typename CSRsp::alloc_type;
  std::vector<CSRsp> B;
  B.reserve(A.size());
  for (auto& v : A)
    B.emplace_back(CSRsp(v, Alloc{v.getAlloc()}));
  return B;
}

} // namespace shm

} // namespace csr

#endif
