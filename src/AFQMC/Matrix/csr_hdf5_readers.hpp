//////////////////////////////////////////////////////////////////////////////////////
//// This file is distributed under the University of Illinois/NCSA Open Source License.
//// See LICENSE file in top directory for details.
////
//// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
////
//// File developed by:
////
//// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_CSR_HDF5_READERS_HPP
#define QMCPLUSPLUS_AFQMC_CSR_HDF5_READERS_HPP


#include <cassert>
#include <complex>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <vector>
#include <numeric>
#if defined(USE_MPI)
#include <mpi.h>
#endif

#include "Configuration.h"
#include "Utilities/FairDivide.h"
#include "hdf/hdf_archive.h"

#include "AFQMC/config.0.h"
#include "AFQMC/Utilities/afqmc_TTI.hpp"
#include "AFQMC/Matrix/array_partition.hpp"
#include "AFQMC/Matrix/matrix_emplace_wrapper.hpp"

#include "mpi3/communicator.hpp"
#include "mpi3/shared_communicator.hpp"
#include "mpi3/shared_window.hpp"

#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Matrix/coo_matrix.hpp"
#include "AFQMC/Utilities/hdf5_consistency_helper.hpp"

#define MORE 1555
#define LAST 2777

namespace qmcplusplus
{
namespace afqmc
{
namespace csr_hdf5
{
using communicator        = boost::mpi3::communicator;
using shared_communicator = boost::mpi3::shared_communicator;

// Careful here!!!
// I disabled copy assignment/constructor in all sparse matrices, so this is meant to rely
// on either an explicit move assignment when called or on RVO!
// Also notice that the communicators can be anything, even self!
template<class SparseArray2D,
         class Alloc,
         typename = typename std::enable_if<std::decay<SparseArray2D>::type::sparse>::type,
         typename = typename std::enable_if<std::decay<SparseArray2D>::type::dimensionality == -2>::type>
inline SparseArray2D HDF2CSR(hdf_archive& dump, shared_communicator& node) //, Alloc alloc_)
{
  using value_type = typename SparseArray2D::value_type;
  using index_type = typename SparseArray2D::index_type;
  using int_type   = typename SparseArray2D::int_type;

  // Need to read:
  // - dims: nrow,ncols, nnz
  // - data_
  // - jdata_
  // - pointers_begin_
  // - pointers_end_

  using size_type = ma::sparse::size_type;

  size_type nrows, ncols, nnz;
  std::vector<size_type> dims(3);
  if (node.root())
  {
    if (!dump.readEntry(dims, "dims"))
      APP_ABORT("Problems reading dims in csr_from_hdf5. \n");
    assert(dims.size() == 3);
  }
  node.broadcast_n(dims.begin(), dims.size());
  nrows = dims[0];
  ncols = dims[1];
  nnz   = dims[2];

  std::vector<int_type> nnz_per_row(nrows);
  std::vector<int_type> ptrb;
  std::vector<int_type> ptre;
  if (node.root())
  {
    if (!dump.readEntry(ptrb, "pointers_begin_"))
      APP_ABORT("Problems reading pointers_begin_ in csr_from_hdf5. \n");
    if (!dump.readEntry(ptre, "pointers_end_"))
      APP_ABORT("Problems reading pointers_end_ in csr_from_hdf5. \n");
    assert(ptrb.size() == nrows);
    assert(ptre.size() == nrows);
    for (size_type i = 0; i < nrows; i++)
      nnz_per_row[i] = ptre[i] - ptrb[i];
  }

  node.broadcast_n(nnz_per_row.begin(), nnz_per_row.size());
  SparseArray2D SpM(std::tuple<std::size_t, std::size_t>{nrows, ncols}, std::tuple<std::size_t, std::size_t>{0, 0},
                    nnz_per_row, Alloc(node));

  if (node.root())
  {
    std::vector<value_type> data;
    std::vector<index_type> jdata;
    if (!dump.readEntry(data, "data_"))
      APP_ABORT("Problems reading data_ in csr_from_hdf5. \n");
    if (!dump.readEntry(jdata, "jdata_"))
      APP_ABORT("Problems reading jdata_ in csr_from_hdf5. \n");
    if (data.size() != nnz)
      APP_ABORT("Problems with data_ array in csr_from_hdf5. \n");
    if (jdata.size() != nnz)
      APP_ABORT("Problems with jdata_ array in csr_from_hdf5. \n");
    size_type cnt = 0;
    for (index_type i = 0; i < static_cast<index_type>(nrows); i++)
    {
      size_type nt = static_cast<size_type>(ptre[i] - ptrb[i]);
      for (size_type nn = 0; nn < nt; ++nn, ++cnt)
        emplace(SpM, std::make_tuple(i, jdata[cnt], data[cnt]));
    }
  }

  node.barrier();

  return SpM;
}

template<class SparseArray2D,
         typename = typename std::enable_if<std::decay<SparseArray2D>::type::sparse>::type,
         typename = typename std::enable_if<std::decay<SparseArray2D>::type::dimensionality == -2>::type>
inline void CSR2HDF(hdf_archive& dump, SparseArray2D const& SpM)
{
  using value_type = typename SparseArray2D::value_type;
  using index_type = typename SparseArray2D::index_type;
  using int_type   = typename SparseArray2D::int_type;

  // Need to write:
  // - dims: nrow,ncols, nnz
  // - data_
  // - jdata_
  // - pointers_begin_
  // - pointers_end_

  using size_type = ma::sparse::size_type;
  size_type nnz   = SpM.num_non_zero_elements();
  size_type nrows = SpM.size(0);
  size_type ncols = SpM.size(1);
  {
    std::vector<size_type> dims{nrows, ncols, nnz};
    dump.write(dims, "dims");
  }
  {
    std::vector<value_type> data(nnz);
    size_type cnt = 0;
    for (size_type i = 0; i < nrows; i++)
    {
      size_type nt = SpM.num_non_zero_elements(i);
      std::copy_n(to_address(SpM.non_zero_values_data(i)), nt, data.data() + cnt);
      cnt += nt;
    }
    dump.write(data, "data_");
  }
  {
    std::vector<index_type> jdata(nnz);
    size_type cnt = 0;
    for (size_type i = 0; i < nrows; i++)
    {
      size_type nt = SpM.num_non_zero_elements(i);
      std::copy_n(to_address(SpM.non_zero_indices2_data(i)), nt, jdata.data() + cnt);
      cnt += nt;
    }
    dump.write(jdata, "jdata_");
  }
  {
    std::vector<int_type> ptrb(nrows);
    std::vector<int_type> ptre(nrows);
    size_type cnt = 0;
    for (size_type i = 0; i < nrows; i++)
    {
      ptrb[i] = cnt;
      cnt += SpM.num_non_zero_elements(i);
      ptre[i] = cnt;
    }
    dump.write(ptrb, "pointers_begin_");
    dump.write(ptre, "pointers_end_");
  }
}

template<class value_type, class index_type, class container, class matrix_map_, class task_group>
inline void multiple_reader_hdf5_csr(container& Q,
                                     matrix_map_ const& map_,
                                     hdf_archive& dump,
                                     task_group& TG,
                                     int n_working_cores)
{
  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID(), coreid = TG.getCoreID();

  // process hdf5 file
  if (coreid < n_working_cores)
  {
    std::vector<int> block_size;
    if (!dump.readEntry(block_size, std::string("block_sizes")))
    {
      app_error() << " Error in multiple_reader_hdf5_csr: Problems reading block_sizes dataset. \n";
      APP_ABORT(" Error in multiple_reader_hdf5_csr: Problems reading block_sizes dataset. \n");
    }

    int nblk  = block_size.size();
    int ntmax = *std::max_element(block_size.begin(), block_size.end());
    std::vector<IndexType> ivec, ivec2;
    std::vector<ValueType> vvec, vvec2;
    ivec.reserve(2 * ntmax);
    vvec.reserve(ntmax);
    ivec2.reserve(2 * ntmax);
    vvec2.reserve(ntmax);

    std::vector<int> pool_dist;
    // divide blocks over n_working_cores groups
    FairDivide(nblk, n_working_cores, pool_dist);

    int first_local_block = pool_dist[coreid];
    int last_local_block  = pool_dist[coreid + 1];
    int nbtot             = last_local_block - first_local_block;
    int niter             = nbtot / nnodes + std::min(nbtot % nnodes, 1);

    for (int iter = 0; iter < niter; iter++)
    {
      int myblock_number = first_local_block + nnodes * iter + nodeid;
      int first_block    = first_local_block + nnodes * iter;
      int last_block     = std::min(first_block + nnodes, last_local_block);
      if (myblock_number < last_local_block)
      {
        ivec.resize(2 * block_size[myblock_number]);
        vvec.resize(block_size[myblock_number]);
        if (block_size[myblock_number] > 0)
        {
          if (!dump.readEntry(ivec, std::string("index_") + std::to_string(myblock_number)))
          {
            app_error() << " Error in multiple_reader_hdf5_csr: Problems reading index_" << myblock_number
                        << " dataset. \n";
            APP_ABORT(" Error in multiple_reader_hdf5_csr: Problems reading index_ dataset. \n");
          }
          if (!readComplexOrReal(dump, std::string("vals_") + std::to_string(myblock_number), vvec))
          {
            app_error() << " Error in multiple_reader_hdf5_csr: Problems reading vals_" << myblock_number
                        << " dataset. \n";
            APP_ABORT(" Error in multiple_reader_hdf5_csr: Problems reading vals_ dataset. \n");
          }
        }
      }
      for (int k = first_block, ipr = 0; k < last_block; k++, ipr++)
      {
        ivec2.resize(2 * block_size[k]);
        vvec2.resize(block_size[k]);
        if (ipr == nodeid)
        {
          assert(myblock_number == k);
          std::copy(ivec.begin(), ivec.end(), ivec2.begin());
          TG.Cores().broadcast_n(ivec2.begin(), ivec2.size(), ipr);
          std::copy(vvec.begin(), vvec.end(), vvec2.begin());
          TG.Cores().broadcast_n(vvec2.begin(), vvec2.size(), ipr);
        }
        else
        {
          TG.Cores().broadcast_n(ivec2.begin(), ivec2.size(), ipr);
          TG.Cores().broadcast_n(vvec2.begin(), vvec2.size(), ipr);
        }

        for (int ik = 0, ikend = block_size[k]; ik < ikend; ik++)
        {
          if (map_(ivec2[2 * ik], ivec2[2 * ik + 1], vvec2[ik]))
          {
            auto i_ = map_.map(ivec2[2 * ik], ivec2[2 * ik + 1]);
            Q.emplace_back(
                std::forward_as_tuple(index_type(i_[0]), index_type(i_[1]), static_cast<value_type>(vvec2[ik])));
          }
        }
      }
    }
  }
}

template<class value_type, class index_type, class container, class matrix_map_, class task_group>
inline void read_csr_matrix_from_hdf_into_distributed_container(container& Q,
                                                                matrix_map_ const& map_,
                                                                hdf_archive& dump,
                                                                task_group& TG)
{
  int rank = TG.Global().rank();
  int size = TG.Global().size();

  std::vector<int> block_size;
  if (!dump.readEntry(block_size, std::string("block_sizes")))
  {
    app_error() << " Error in multiple_reader_hdf5_csr: Problems reading block_sizes dataset. \n";
    APP_ABORT(" Error in multiple_reader_hdf5_csr: Problems reading block_sizes dataset. \n");
  }

  int nblk = block_size.size();
  int first_block, last_block;
  std::tie(first_block, last_block) = FairDivideBoundary(rank, nblk, size);

  if (first_block >= last_block)
    return;

  int ntmax = *std::max_element(block_size.begin() + first_block, block_size.begin() + last_block);
  std::vector<IndexType> ivec;
  std::vector<ValueType> vvec;
  ivec.reserve(2 * ntmax);
  vvec.reserve(ntmax);

  for (int k = first_block, ipr = 0; k < last_block; k++, ipr++)
  {
    ivec.resize(2 * block_size[k]);
    vvec.resize(block_size[k]);
    if (block_size[k] > 0)
    {
      if (!dump.readEntry(ivec, std::string("index_") + std::to_string(k)))
      {
        app_error() << " Error in multiple_reader_hdf5_csr: Problems reading index_" << k << " dataset. \n";
        APP_ABORT(" Error in multiple_reader_hdf5_csr: Problems reading index_ dataset. \n");
      }
      if (!readComplexOrReal(dump, std::string("vals_") + std::to_string(k), vvec))
      {
        app_error() << " Error in multiple_reader_hdf5_csr: Problems reading vals_" << k << " dataset. \n";
        APP_ABORT(" Error in multiple_reader_hdf5_csr: Problems reading vals_ dataset. \n");
      }
      for (int ik = 0, ikend = block_size[k]; ik < ikend; ik++)
      {
        if (map_(ivec[2 * ik], ivec[2 * ik + 1], vvec[ik]))
        {
          auto i_ = map_.map(ivec[2 * ik], ivec[2 * ik + 1]);
          Q.emplace_back(
              std::forward_as_tuple(index_type(i_[0]), index_type(i_[1]), static_cast<value_type>(vvec[ik])));
        }
      }
    }
  }
}

/* 
 * In local_count, split is assumed to be the same on all cores of a node only, 
 *       so blocks on different nodes are communicated and accumulators are reduced over all cores 
 *       on a node. 
 */
template<class matrix_partition, class IVec, class task_group>
inline void multiple_reader_local_count(hdf_archive& dump,
                                        matrix_partition const& split,
                                        IVec& counts,
                                        task_group& TG,
                                        int n_working_cores)
{
  assert(counts.size() == split.range());
  std::fill(counts.begin(), counts.end(), 0);
  assert(sizeof(counts[0]) == sizeof(int));

  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID(), coreid = TG.getCoreID();

  // process hdf5 file
  if (coreid < n_working_cores)
  {
    std::vector<int> block_size;
    if (!dump.readEntry(block_size, std::string("block_sizes")))
    {
      app_error() << " Error in multiple_reader_count_entries: Problems reading block_sizes dataset. \n";
      APP_ABORT(" Error in multiple_reader_count_entries: Problems reading block_sizes dataset. \n");
    }

    int ntmax = *std::max_element(block_size.begin(), block_size.end());
    std::vector<IndexType> ivec, ivec2;
    std::vector<ValueType> vvec, vvec2;
    ivec.reserve(2 * ntmax);
    vvec.reserve(ntmax);
    ivec2.reserve(2 * ntmax);
    vvec2.reserve(ntmax);

    int nblk = block_size.size();
    std::vector<int> pool_dist;
    // divide blocks over n_working_cores groups
    FairDivide(nblk, n_working_cores, pool_dist);

    int first_local_block = pool_dist[coreid];
    int last_local_block  = pool_dist[coreid + 1];
    int nbtot             = last_local_block - first_local_block;
    int niter             = nbtot / nnodes + std::min(nbtot % nnodes, 1);

    for (int iter = 0; iter < niter; iter++)
    {
      int myblock_number = first_local_block + nnodes * iter + nodeid;
      int first_block    = first_local_block + nnodes * iter;
      int last_block     = std::min(first_block + nnodes, last_local_block);
      if (myblock_number < last_local_block)
      {
        ivec.resize(2 * block_size[myblock_number]);
        vvec.resize(block_size[myblock_number]);
        if (block_size[myblock_number] > 0)
        {
          if (!dump.readEntry(ivec, std::string("index_") + std::to_string(myblock_number)))
          {
            app_error() << " Error in multiple_reader_hdf5_SpMat: Problems reading index_" << myblock_number
                        << " dataset. \n";
            APP_ABORT(" Error in multiple_reader_hdf5_SpMat: Problems reading index dataset. \n");
          }
          if (!readComplexOrReal(dump, std::string("vals_") + std::to_string(myblock_number), vvec))
          {
            app_error() << " Error in multiple_reader_hdf5_SpMat: Problems reading vals_" << myblock_number
                        << " dataset. \n";
            APP_ABORT(" Error in multiple_reader_hdf5_SpMat: Problems reading vals_ dataset. \n");
          }
        }
      }
      for (int k = first_block, ipr = 0; k < last_block; k++, ipr++)
      {
        ivec2.resize(2 * block_size[k]);
        vvec2.resize(block_size[k]);
        if (ipr == nodeid)
        {
          assert(myblock_number == k);
          std::copy(ivec.begin(), ivec.end(), ivec2.begin());
          TG.Cores().broadcast_n(ivec2.begin(), ivec2.size(), ipr);
          std::copy(vvec.begin(), vvec.end(), vvec2.begin());
          TG.Cores().broadcast_n(vvec2.begin(), vvec2.size(), ipr);
        }
        else
        {
          TG.Cores().broadcast_n(ivec2.begin(), ivec2.size(), ipr);
          TG.Cores().broadcast_n(vvec2.begin(), vvec2.size(), ipr);
        }
        for (int ik = 0, ikend = block_size[k]; ik < ikend; ik++)
        {
          auto ip = split.map(ivec2[2 * ik], ivec2[2 * ik + 1], vvec2[ik]);
          if (ip >= 0)
            counts[ip]++;
        }
      } // k < last_block
    }
  }

  TG.Node().all_reduce_in_place_n(counts.begin(), counts.size(), std::plus<>());
}

/*
 * In global_count, all cores accumulate entries the same quantity. 
 *       Split is expected to be the same on all nodes and no communication of blocks is needed or 
 *       performed. Final accumulators are reduced over the global communicator.
 */
template<class matrix_partition, class IVec, class task_group>
inline void multiple_reader_global_count(hdf_archive& dump,
                                         matrix_partition const& split,
                                         IVec& counts,
                                         task_group& TG,
                                         int n_working_cores)
{
  assert(counts.size() == split.range());
  std::fill(counts.begin(), counts.end(), 0);
  assert(sizeof(counts[0]) == sizeof(int));

  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID(), coreid = TG.getCoreID();

  // process hdf5 file
  if (coreid < n_working_cores)
  {
    std::vector<int> block_size;
    if (!dump.readEntry(block_size, std::string("block_sizes")))
    {
      app_error() << " Error in multiple_reader_count_entries: Problems reading ***_block_sizes dataset. \n";
      APP_ABORT(" Error in multiple_reader_count_entries: Problems reading ***_block_sizes dataset. \n");
    }

    int ntmax = *std::max_element(block_size.begin(), block_size.end());
    std::vector<IndexType> ivec;
    std::vector<ValueType> vvec;
    ivec.reserve(2 * ntmax);
    vvec.reserve(ntmax);

    int nblk = block_size.size();
    std::vector<int> pool_dist;
    // divide blocks over n_working_cores groups
    FairDivide(nblk, n_working_cores, pool_dist);

    int first_local_block = pool_dist[coreid];
    int last_local_block  = pool_dist[coreid + 1];

    for (int ib = first_local_block; ib < last_local_block; ib++)
    {
      if (ib % nnodes != nodeid)
        continue;
      ivec.resize(2 * block_size[ib]);
      vvec.resize(block_size[ib]);
      if (block_size[ib] > 0)
      {
        if (!dump.readEntry(ivec, std::string("index_") + std::to_string(ib)))
        {
          app_error() << " Error in multiple_reader_count_entries: Problems reading index_" << ib << " dataset. \n";
          APP_ABORT(" Error in multiple_reader_count_entries: Problems reading index_  dataset. \n");
        }
        if (!readComplexOrReal(dump, std::string("vals_") + std::to_string(ib), vvec))
        {
          app_error() << " Error in multiple_reader_count_entries: Problems reading vals_" << ib << " dataset. \n";
          APP_ABORT(" Error in multiple_reader_count_entries: Problems reading vals_ dataset. \n");
        }
        for (int ik = 0, ikend = block_size[ib]; ik < ikend; ik++)
        {
          int ip = split.map(ivec[2 * ik], ivec[2 * ik + 1], vvec[ik]);
          if (ip > 0)
            assert(ip < counts.size());
          if (ip >= 0)
            counts[ip]++;
        }
      }
    }
  }

  TG.Global().all_reduce_in_place_n(counts.begin(), counts.size(), std::plus<>());
}

// right now exclusive for shm case
template<class SparseArray2D,
         class task_group_,
         typename = typename std::enable_if<std::decay<SparseArray2D>::type::sparse>::type,
         typename = typename std::enable_if<std::decay<SparseArray2D>::type::dimensionality == -2>::type>
inline SparseArray2D column_distributed_CSR_from_HDF(hdf_archive& dump,
                                                     task_group_& TG,
                                                     int nrows,
                                                     int ncols,
                                                     double cutoff1bar,
                                                     int nread = -1)
{
  using value_type = typename SparseArray2D::value_type;
  using index_type = typename SparseArray2D::index_type;
  using int_type   = typename SparseArray2D::int_type;

  // taken from Hamiltonians/HamiltonianFactory_Helper::read_V2fact
  using counter = qmcplusplus::afqmc::sparse_matrix_element_counter;
  using mat_map = qmcplusplus::afqmc::matrix_map;
  // Take Alloc and is_root from SparseArray2D itself
  using Alloc = shared_allocator<value_type>;
  using ucsr_matrix =
      ma::sparse::ucsr_matrix<value_type, index_type, int_type, shared_allocator<value_type>, ma::sparse::is_root>;

  // if not specified, every core reads
  if (nread <= 0)
    nread = TG.Node().size();

  int c0, cN;
  bool distribute_Ham = (TG.getNGroupsPerTG() > 1);
  std::vector<int> row_counts(nrows);

  // calculate column range that belong to this node
  if (distribute_Ham)
  {
    // count number of non-zero elements along each column (Chol Vec)
    std::vector<IndexType> col_counts(ncols);
    std::vector<IndexType> sets(TG.getNGroupsPerTG() + 1);
    csr_hdf5::multiple_reader_global_count(dump, counter(false, nrows, ncols, 0, nrows, 0, ncols, cutoff1bar),
                                           col_counts, TG, nread);

    if (TG.Global().root())
    {
      simple_matrix_partition<task_group_, IndexType, RealType> split(nrows, ncols, cutoff1bar);
      split.partition(TG, false, col_counts, sets);

      app_log() << " Partitioning of columns in column_distributed_CSR_from_HDF: ";
      for (int i = 0; i <= TG.getNGroupsPerTG(); i++)
        app_log() << sets[i] << " ";
      app_log() << std::endl;
      app_log() << " Number of terms in each partitioning: ";
      for (int i = 0; i < TG.getNGroupsPerTG(); i++)
        app_log() << accumulate(col_counts.begin() + sets[i], col_counts.begin() + sets[i + 1], 0) << " ";
      app_log() << std::endl;
    }

    TG.Global().broadcast(sets.begin(), sets.end());
    c0 = sets[TG.getLocalGroupNumber()];
    cN = sets[TG.getLocalGroupNumber() + 1];

    csr_hdf5::multiple_reader_local_count(dump, counter(true, nrows, ncols, 0, nrows, c0, cN, cutoff1bar), row_counts,
                                          TG, nread);
  }
  else
  {
    c0 = 0;
    cN = ncols;

    // should be faster if ham is not distributed
    csr_hdf5::multiple_reader_global_count(dump, counter(true, nrows, ncols, 0, nrows, 0, ncols, cutoff1bar),
                                           row_counts, TG, nread);
  }

  ucsr_matrix ucsr(tp_ul_ul{nrows, cN - c0}, tp_ul_ul{0, c0}, row_counts, Alloc(TG.Node()));
  csr::matrix_emplace_wrapper<ucsr_matrix> csr_wrapper(ucsr, TG.Node());

  csr_hdf5::multiple_reader_hdf5_csr<value_type, index_type>(csr_wrapper,
                                                             mat_map(false, true, nrows, ncols, 0, nrows, c0, cN,
                                                                     cutoff1bar),
                                                             dump, TG, nread);
  csr_wrapper.push_buffer();
  TG.node_barrier();

  return SparseArray2D(std::move(ucsr));
}

// right now exclusive for shm case
template<class SparseArray2D,
         class task_group_,
         typename = typename std::enable_if<std::decay<SparseArray2D>::type::sparse>::type,
         typename = typename std::enable_if<std::decay<SparseArray2D>::type::dimensionality == -2>::type>
inline SparseArray2D unstructured_distributed_CSR_from_HDF(hdf_archive& dump,
                                                           task_group_& TG,
                                                           int nrows,
                                                           int ncols,
                                                           double cutoff1bar,
                                                           int nread = -1)
{
  using value_type = typename SparseArray2D::value_type;
  using index_type = typename SparseArray2D::index_type;
  using int_type   = typename SparseArray2D::int_type;

  using mat_map = qmcplusplus::afqmc::matrix_map;
  // Take Alloc and is_root from SparseArray2D itself
  using Alloc = shared_allocator<value_type>;
  using ucsr_matrix =
      ma::sparse::ucsr_matrix<value_type, index_type, int_type, shared_allocator<value_type>, ma::sparse::is_root>;

  bool distribute_Ham = (TG.getNGroupsPerTG() > 1);

  if (distribute_Ham)
  {
    // if not specified, every core reads
    if (nread <= 0)
      nread = TG.Node().size();

    using tvec = std::vector<std::tuple<index_type, index_type, value_type>>;
    tvec Q;
    if (TG.Node().rank() < nread)
      Q.reserve(100000); // reserve some space

    csr_hdf5::read_csr_matrix_from_hdf_into_distributed_container<value_type, index_type>(Q,
                                                                                          mat_map(false, false, nrows,
                                                                                                  ncols, 0, nrows, 0,
                                                                                                  ncols, cutoff1bar),
                                                                                          dump, TG);

    return csr::shm::construct_distributed_csr_matrix_from_distributed_containers<SparseArray2D>(Q, nrows, ncols, TG);
  }
  else
  {
    return column_distributed_CSR_from_HDF<SparseArray2D>(dump, TG, nrows, ncols, cutoff1bar, nread);
  }
}

template<class SparseArray2D,
         class task_group_,
         typename = typename std::enable_if<std::decay<SparseArray2D>::type::sparse>::type,
         typename = typename std::enable_if<std::decay<SparseArray2D>::type::dimensionality == -2>::type>
inline void write_distributed_CSR_to_HDF(SparseArray2D const& SpM, hdf_archive& dump, task_group_& TG)
{
  using value_type = typename SparseArray2D::value_type;
  using index_type = typename SparseArray2D::index_type;
  using int_type   = typename SparseArray2D::int_type;
  using std::size_t;

  size_t nnz = SpM.num_non_zero_elements();
  if (nnz == size_t(0))
    return;

  if (TG.getNGroupsPerTG() > 1)
  {
    if (TG.Global().root())
    {
      std::vector<int> block_sizes;
      block_sizes.reserve((nnz - size_t(1)) / CSR_HDF_BLOCK_SIZE + 1);
      std::vector<value_type> vvec;
      vvec.reserve(CSR_HDF_BLOCK_SIZE);
      std::vector<index_type> ivec;
      ivec.reserve(2 * CSR_HDF_BLOCK_SIZE);
      size_t cnt = 0;
      for (size_t r = 0; r < SpM.size(0); ++r)
      {
        size_t nzr = SpM.num_non_zero_elements(r);
        auto val   = SpM.non_zero_values_data(r);
        auto col   = SpM.non_zero_indices2_data(r);
        for (size_t c = 0; c < nzr; ++c)
        {
          vvec.push_back(val[c]);
          ivec.push_back(index_type(r));
          ivec.push_back(col[c]);
          if (vvec.size() == CSR_HDF_BLOCK_SIZE)
          {
            dump.write(vvec, std::string("vals_") + std::to_string(cnt));
            dump.write(ivec, std::string("index_") + std::to_string(cnt++));
            ivec.clear();
            vvec.clear();
            block_sizes.push_back(CSR_HDF_BLOCK_SIZE);
          }
        }
      }
      if (vvec.size() > 0)
      {
        block_sizes.push_back(vvec.size());
        dump.write(vvec, std::string("vals_") + std::to_string(cnt));
        dump.write(ivec, std::string("index_") + std::to_string(cnt++));
        ivec.clear();
        vvec.clear();
      }

      int nnodes_per_TG = TG.getNGroupsPerTG();
      for (int core = 1; core < nnodes_per_TG; ++core)
      {
        size_t nnz_;
        TG.Cores().receive_n(&nnz_, 1, core, core);
        int nblk = ((nnz_ - size_t(1)) / CSR_HDF_BLOCK_SIZE + 1);
        for (int ni = 0; ni < nblk; ++ni)
        {
          size_t sz = std::min(size_t(CSR_HDF_BLOCK_SIZE), nnz_ - size_t(ni * CSR_HDF_BLOCK_SIZE));
          vvec.resize(sz);
          ivec.resize(2 * sz);
          TG.Cores().receive_n(vvec.data(), sz, core, nnodes_per_TG + core);
          TG.Cores().receive_n(ivec.data(), 2 * sz, core, 2 * nnodes_per_TG + core);
          dump.write(vvec, std::string("vals_") + std::to_string(cnt));
          dump.write(ivec, std::string("index_") + std::to_string(cnt++));
          ivec.clear();
          vvec.clear();
          block_sizes.push_back(sz);
        }
      }

      dump.write(block_sizes, "block_sizes");
    }
    else if (TG.Node().root() && TG.Cores().rank() < TG.getNGroupsPerTG())
    { // only one TG writes
      int nnodes_per_TG = TG.getNGroupsPerTG();
      int rank          = TG.Cores().rank();
      TG.Cores().send_n(&nnz, 1, 0, rank);
      std::vector<value_type> vvec;
      vvec.reserve(CSR_HDF_BLOCK_SIZE);
      std::vector<index_type> ivec;
      ivec.reserve(2 * CSR_HDF_BLOCK_SIZE);
      index_type offset_r = index_type(SpM.global_origin()[0]);
      index_type offset_c = index_type(SpM.global_origin()[1]);
      for (size_t r = 0; r < SpM.size(0); ++r)
      {
        size_t nzr = SpM.num_non_zero_elements(r);
        auto val   = SpM.non_zero_values_data(r);
        auto col   = SpM.non_zero_indices2_data(r);
        for (size_t c = 0; c < nzr; ++c)
        {
          vvec.push_back(val[c]);
          ivec.push_back(index_type(r) + offset_r);
          ivec.push_back(col[c] + offset_c);
          if (vvec.size() == CSR_HDF_BLOCK_SIZE)
          {
            TG.Cores().send_n(vvec.data(), vvec.size(), 0, nnodes_per_TG + rank);
            TG.Cores().send_n(ivec.data(), ivec.size(), 0, 2 * nnodes_per_TG + rank);
            ivec.clear();
            vvec.clear();
          }
        }
      }
      if (vvec.size() > 0)
      {
        TG.Cores().send_n(vvec.data(), vvec.size(), 0, nnodes_per_TG + rank);
        TG.Cores().send_n(ivec.data(), ivec.size(), 0, 2 * nnodes_per_TG + rank);
      }
    }
  }
  else
  {
    if (TG.Global().root())
    {
      std::vector<int> block_sizes;
      block_sizes.reserve((nnz - size_t(1)) / CSR_HDF_BLOCK_SIZE + 1);
      std::vector<value_type> vvec;
      vvec.reserve(CSR_HDF_BLOCK_SIZE);
      std::vector<index_type> ivec;
      ivec.reserve(2 * CSR_HDF_BLOCK_SIZE);
      size_t cnt = 0;
      for (size_t r = 0; r < SpM.size(0); ++r)
      {
        size_t nzr = SpM.num_non_zero_elements(r);
        auto val   = SpM.non_zero_values_data(r);
        auto col   = SpM.non_zero_indices2_data(r);
        for (size_t c = 0; c < nzr; ++c)
        {
          vvec.push_back(val[c]);
          ivec.push_back(index_type(r));
          ivec.push_back(col[c]);
          if (vvec.size() == CSR_HDF_BLOCK_SIZE)
          {
            dump.write(vvec, std::string("vals_") + std::to_string(cnt));
            dump.write(ivec, std::string("index_") + std::to_string(cnt++));
            ivec.clear();
            vvec.clear();
            block_sizes.push_back(CSR_HDF_BLOCK_SIZE);
          }
        }
      }
      if (vvec.size() > 0)
      {
        block_sizes.push_back(vvec.size());
        dump.write(vvec, std::string("vals_") + std::to_string(cnt));
        dump.write(ivec, std::string("index_") + std::to_string(cnt++));
      }
      assert(cnt == ((nnz - size_t(1)) / CSR_HDF_BLOCK_SIZE + 1));
      dump.write(block_sizes, "block_sizes");
    }
  }
  TG.Global().barrier();
}

} // namespace csr_hdf5

} // namespace afqmc

} // namespace qmcplusplus

#endif
