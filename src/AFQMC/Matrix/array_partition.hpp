//////////////////////////////////////////////////////////////////////////////////////
////// This file is distributed under the University of Illinois/NCSA Open Source License.
////// See LICENSE file in top directory for details.
//////
////// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//////
////// File developed by:
//////
////// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_ARRAY_PARTITION_HPP
#define QMCPLUSPLUS_AFQMC_ARRAY_PARTITION_HPP

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
#include "Message/CommOperators.h"

#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/config.0.h"

namespace qmcplusplus
{
namespace afqmc
{
/*
 * Object used to count number of non-zero elements in a stored (e.g. hdf5) sparse matrix.
 * Given a subset of the global matrix, the class implements a map from {i,j}->k,
 * where k is an integer between 0:counter_range used to count elements.
 * Either row or column reductions are implemented.
 */
struct sparse_matrix_element_counter
{
  using IType = std::size_t;
  using DType = double;

public:
  sparse_matrix_element_counter(bool row_,
                                IType ngr_,
                                IType ngc_,
                                IType r0_,
                                IType r1_,
                                IType c0_,
                                IType c1_,
                                DType cut_)
      : byrow_(row_), ngr(ngr_), ngc(ngc_), r0(r0_), r1(r1_), c0(c0_), c1(c1_), cut(std::abs(cut_))
  {
    assert(ngr > 0);
    assert(ngc > 0);
    assert(r0 >= 0);
    assert(c0 >= 0);
    assert(r1 > r0);
    assert(c1 > c0);
    assert(r1 <= ngr);
    assert(c1 <= ngc);
  }

  auto range() const { return ((byrow_) ? (r1 - r0) : (c1 - c0)); }

  template<class integer_, class double_>
  int map(integer_ i, integer_ j, double_ v) const
  {
    if (std::abs(v) <= cut || IType(i) < r0 || IType(i) >= r1 || IType(j) < c0 || IType(j) >= c1)
      return -1;
    else
      return (byrow_) ? (IType(i) - r0) : (IType(j) - c0);
  }

private:
  bool byrow_;
  IType ngr, ngc, r0, r1, c0, c1;
  DType cut;
};

struct matrix_map
{
  using IType = std::size_t;
  using DType = double;

public:
  matrix_map(bool row_, bool col_, IType ngr_, IType ngc_, IType r0_, IType r1_, IType c0_, IType c1_, DType cut_)
      : shiftrow_(row_), shiftcol_(col_), ngr(ngr_), ngc(ngc_), r0(r0_), r1(r1_), c0(c0_), c1(c1_), cut(std::abs(cut_))
  {
    assert(ngr > 0);
    assert(ngc > 0);
    assert(r0 >= 0);
    assert(c0 >= 0);
    assert(r1 > r0);
    assert(c1 > c0);
    assert(r1 <= ngr);
    assert(c1 <= ngc);
  }

  auto range() const { return std::array<IType, 2>{r1 - r0, c1 - c0}; }

  template<class integer_, class double_ = double>
  bool operator()(integer_ i, integer_ j, double_ v = double_(0.0)) const
  {
    if (std::abs(v) <= cut || IType(i) < r0 || IType(i) >= r1 || IType(j) < c0 || IType(j) >= c1)
      return false;
    return true;
  }

  template<class integer_>
  std::array<integer_, 2> map(integer_ i, integer_ j) const
  {
    return {(shiftrow_ ? integer_(IType(i) - r0) : i), (shiftcol_ ? integer_(IType(j) - c0) : j)};
  }

private:
  bool shiftrow_, shiftcol_;
  IType ngr, ngc, r0, r1, c0, c1;
  DType cut;
};

/*
 * object that encapsulates the partitioning of a 2-D array (matrix) 
 */
template<class task_group, typename IType, typename DType>
struct simple_matrix_partition
{
  simple_matrix_partition(IType nr, IType nc, DType rc = 0)
      : cut(std::abs(rc)), r0(0), r1(nr), c0(0), c1(nc), nrows(nr), ncols(nc)
  {}

  ~simple_matrix_partition() {}

  template<typename DType2>
  inline bool local(const IType i, const IType j, const DType2 v) const
  {
    return (i >= r0) && (i < r1) && (j >= c0) && (j < c1) && (std::abs(v) > cut);
  }

  // this breaks a matrix dimension over TG.nnodes, assuming homogeneous blocks
  inline void partition(const task_group& TG, bool byRow, std::vector<IType>& sets)
  {
    int nnodes = TG.getNGroupsPerTG();
    IType cnt  = 0;
    int nblk;
    if (byRow)
      nblk = nrows;
    else
      nblk = ncols;
    assert(nblk >= nnodes);
    sets.resize(nnodes + 1);
    sets[0] = 0;
    if (nnodes > 1)
    {
      FairDivide(nblk, nnodes, sets);
      int node_number = TG.getLocalGroupNumber();
      if (byRow)
      {
        r0 = sets[node_number];
        r1 = sets[node_number + 1];
      }
      else
      {
        c0 = sets[node_number];
        c1 = sets[node_number + 1];
      }
    }
    else
    {
      if (byRow)
      {
        r0      = 0;
        r1      = nrows;
        sets[1] = nrows;
      }
      else
      {
        c0      = 0;
        c1      = ncols;
        sets[1] = ncols;
      }
    }
  }

  // this breaks a matrix dimension over TG.nnodes
  inline void partition(const task_group& TG, bool byRow, const std::vector<IType>& counts, std::vector<IType>& sets)
  {
    int nnodes = TG.getNGroupsPerTG();
    int nblk   = counts.size();
    IType cnt  = 0;
    assert(nblk >= nnodes);
    if (byRow)
      assert(nblk == nrows);
    else
      assert(nblk == ncols);
    sets.resize(nnodes + 1);
    sets[0] = 0;
    if (nnodes > 1)
    {
      std::vector<IType> nv(counts.size() + 1);
      nv[0]                                     = 0;
      typename std::vector<IType>::iterator itn = nv.begin() + 1;
      for (typename std::vector<IType>::const_iterator itc = counts.begin(), ite = counts.end(); itc != ite;
           itc++, itn++)
      {
        cnt += *(itc);
        (*itn) = cnt;
      }
      balance_partition_ordered_set(counts.size(), nv.data(), sets);
      int node_number = TG.getLocalGroupNumber();
      if (byRow)
      {
        r0 = sets[node_number];
        r1 = sets[node_number + 1];
      }
      else
      {
        c0 = sets[node_number];
        c1 = sets[node_number + 1];
      }
    }
    else
    {
      if (byRow)
      {
        r0      = 0;
        r1      = nrows;
        sets[1] = nrows;
      }
      else
      {
        c0      = 0;
        c1      = ncols;
        sets[1] = ncols;
      }
    }
  }

  // this breaks a matrix dimension over TG.
  inline void partition_over_TGs(const task_group& TG,
                                 bool byRow,
                                 const std::vector<IType>& counts,
                                 std::vector<IType>& sets)
  {
    int ngrps = TG.getNumberOfTGs();
    int nblk  = counts.size();
    IType cnt = 0;
    assert(nblk >= ngrps);
    if (byRow)
      assert(nblk == nrows);
    else
      assert(nblk == ncols);
    sets.resize(ngrps + 1);
    sets[0] = 0;
    if (ngrps > 1)
    {
      std::vector<IType> nv(counts.size() + 1);
      nv[0]                                     = 0;
      typename std::vector<IType>::iterator itn = nv.begin() + 1;
      for (typename std::vector<IType>::const_iterator itc = counts.begin(), ite = counts.end(); itc != ite;
           itc++, itn++)
      {
        cnt += *(itc);
        (*itn) = cnt;
      }
      balance_partition_ordered_set(counts.size(), nv.data(), sets);
      int node_number = TG.getTGNumber();
      if (byRow)
      {
        r0 = sets[node_number];
        r1 = sets[node_number + 1];
      }
      else
      {
        c0 = sets[node_number];
        c1 = sets[node_number + 1];
      }
    }
    else
    {
      if (byRow)
      {
        r0      = 0;
        r1      = nrows;
        sets[1] = nrows;
      }
      else
      {
        c0      = 0;
        c1      = ncols;
        sets[1] = ncols;
      }
    }
  }

  // this breaks a local segment over TG.ncores_per_TG, assumes homogeneous blocks
  inline void sub_partition(const task_group& TG, bool byRow, std::vector<IType>& sets)
  {
    int ncores = TG.getNCoresPerTG();
    int nblk;
    if (byRow)
      nblk = r1 - r0;
    else
      nblk = c1 - c0;
    assert(nblk >= ncores);
    sets.resize(ncores + 1);
    sets[0] = 0;
    if (ncores > 1)
    {
      FairDivide(nblk, ncores, sets);
      int core_rank = TG.getCoreRank();
      if (byRow)
      {
        sr0 = r0 + sets[core_rank];
        sr1 = r0 + sets[core_rank + 1];
      }
      else
      {
        sc0 = c0 + sets[core_rank];
        sc1 = c0 + sets[core_rank + 1];
      }
    }
    else
    {
      if (byRow)
      {
        sr0     = r0;
        sr1     = r1;
        sets[1] = r1;
      }
      else
      {
        sc0     = c0;
        sc1     = c1;
        sets[1] = c1;
      }
    }
  }

  // this breaks a local segment over TG.ncores_per_TG
  inline void sub_partition(const task_group& TG,
                            bool byRow,
                            const std::vector<IType>& counts,
                            std::vector<IType>& sets)
  {
    int ncores = TG.getNCoresPerTG();
    int nblk   = counts.size();
    IType cnt  = 0;
    assert(nblk >= ncores);
    if (byRow)
      assert(nblk == r1 - r0);
    else
      assert(nblk == c1 - c0);
    sets.resize(ncores + 1);
    sets[0] = 0;
    if (ncores > 1)
    {
      std::vector<IType> nv(counts.size() + 1);
      nv[0]                                     = 0;
      typename std::vector<IType>::iterator itn = nv.begin() + 1;
      for (typename std::vector<IType>::const_iterator itc = counts.begin(), ite = counts.end(); itc != ite;
           itc++, itn++)
      {
        cnt += *(itc);
        (*itn) = cnt;
      }
      balance_partition_ordered_set(counts.size(), nv.data(), sets);
      int core_rank = TG.getCoreRank();
      if (byRow)
      {
        sr0 = r0 + sets[core_rank];
        sr1 = r0 + sets[core_rank + 1];
      }
      else
      {
        sc0 = c0 + sets[core_rank];
        sc1 = c0 + sets[core_rank + 1];
      }
    }
    else
    {
      if (byRow)
      {
        sr0     = r0;
        sr1     = r1;
        sets[1] = r1;
      }
      else
      {
        sc0     = c0;
        sc1     = c1;
        sets[1] = c1;
      }
    }
  }

  std::tuple<int, int, int, int> getLocalPartition() const { return std::make_tuple(r0, r1, c0, c1); }

  std::tuple<int, int, int, int> getSubPartition() const { return std::make_tuple(sr0, sr1, sc0, sc1); }

  DType getCutoff() const { return cut; }

  std::tuple<int, int> getDimensions() const { return std::make_tuple(nrows, ncols); }

private:
  DType cut;
  IType r0, r1;   // lower and upper bond of row segment
  IType c0, c1;   // lower and upper bound of column segment
  IType sr0, sr1; // lower and upper bond of row sub-partition (relative to 0)
  IType sc0, sc1; // lower and upper bound of column sub-partition (relative to 0)
  IType nrows, ncols;
};


} // namespace afqmc

} // namespace qmcplusplus

#endif
