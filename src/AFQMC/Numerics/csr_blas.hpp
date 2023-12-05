//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef CSR_BLAS_HPP
#define CSR_BLAS_HPP

#include <utility> //std::enable_if
#include <cassert>
#include <iostream>
#include <tuple>
#include <complex>
#include "Utilities/FairDivide.h"
#include "AFQMC/Numerics/detail/utilities.hpp"
#include "AFQMC/Matrix/csr_matrix.hpp"

#include <type_traits> // enable_if

using std::complex;
using std::tuple;

namespace csr
{
// y = a*x + y
// y satisfies a 1D multi-array concept
// x satisfies a 1D sparse-array concept
template<class T,
         class SparseArray1D,
         class MultiArray1D,
         typename = typename std::enable_if<std::decay<SparseArray1D>::type::sparse>::type,
         typename = typename std::enable_if<std::decay<SparseArray1D>::type::dimensionality == 1>::type,
         typename = typename std::enable_if<std::decay<MultiArray1D>::type::dimensionality == 1>::type>
MultiArray1D axpy(char TA, T a, SparseArray1D&& x, MultiArray1D&& y)
{
  using ma::conj;
  assert(std::get<0>(x.sizes()) == std::get<0>(y.sizes()));
  auto vals = x.non_zero_values_data();
  auto cols = x.non_zero_indices2_data();
  if (TA == 'C')
    for (std::size_t i = 0, iend = x.num_non_zero_elements(); i < iend; ++i, ++vals, ++cols)
    {
      assert(*cols >= 0 && *cols < y.size());
      y[*cols] += ma::conj(*vals) * a;
    }
  else
    for (std::size_t i = 0, iend = x.num_non_zero_elements(); i < iend; ++i, ++vals, ++cols)
    {
      assert(*cols >= 0 && *cols < y.size());
      y[*cols] += (*vals) * a;
    }
  return std::forward<MultiArray1D>(y);
}

// y = a*x + b*y
// y satisfies a 1D multi-array concept
// x satisfies a 1D sparse-array concept
template<class T,
         class SparseArray1D,
         class MultiArray1D,
         typename = typename std::enable_if<std::decay<SparseArray1D>::type::sparse>::type,
         typename = typename std::enable_if<std::decay<SparseArray1D>::type::dimensionality == 1>::type,
         typename = typename std::enable_if<std::decay<MultiArray1D>::type::dimensionality == 1>::type>
MultiArray1D axpby(char TA, T a, SparseArray1D&& x, T b, MultiArray1D&& y)
{
  using ma::conj;
  assert(x.size(0) == y.size(0));
  auto vals = x.non_zero_values_data();
  auto cols = x.non_zero_indices2_data();
  for (auto& yi : y)
    yi *= b;
  if (TA == 'C')
    for (std::size_t i = 0, iend = x.num_non_zero_elements(); i < iend; ++i, ++vals, ++cols)
    {
      assert(*cols >= 0 && *cols < y.size(0));
      y[*cols] += ma::conj(*vals) * a;
    }
  else
    for (std::size_t i = 0, iend = x.num_non_zero_elements(); i < iend; ++i, ++vals, ++cols)
    {
      assert(*cols >= 0 && *cols < y.size(0));
      y[*cols] += (*vals) * a;
    }
  return std::forward<MultiArray1D>(y);
}

// Dot product between 2 sparse vectors
template<class T, class integer, class VPtr, class JPtr>
inline T csrvv(char TA, char TB, std::tuple<integer, VPtr, JPtr> const& V1, std::tuple<integer, VPtr, JPtr> const& V2)
{
  assert((TA == 'N' or TA == 'C'));
  assert((TB == 'N' or TB == 'C'));
  using std::get;
  T res     = T(0);
  integer i = 0, j = 0;
  const integer n1 = get<0>(V1);
  const VPtr A1    = get<1>(V1);
  const JPtr indx1 = get<2>(V1);
  const integer n2 = get<0>(V2);
  const VPtr A2    = get<1>(V2);
  const JPtr indx2 = get<2>(V2);

  if (TA == 'N' && TB == 'N')
  {
    while (i < n1 && j < n2)
    {
      if (*(indx1 + i) < *(indx2 + j))
        ++i;
      else if (*(indx2 + j) < *(indx1 + i))
        ++j;
      else
      {
        res += *(A1 + i) * (*(A2 + j));
        ++i;
        ++j;
      }
    }
  }
  else if (TA == 'N' && TB == 'C')
  {
    while (i < n1 && j < n2)
    {
      if (*(indx1 + i) < *(indx2 + j))
        ++i;
      else if (*(indx2 + j) < *(indx1 + i))
        ++j;
      else
      {
        res += *(A1 + i) * ma::conj(*(A2 + j));
        ++i;
        ++j;
      }
    }
  }
  else if (TA == 'C' && TB == 'N')
  {
    while (i < n1 && j < n2)
    {
      if (*(indx1 + i) < *(indx2 + j))
        ++i;
      else if (*(indx2 + j) < *(indx1 + i))
        ++j;
      else
      {
        res += ma::conj(*(A1 + i)) * (*(A2 + j));
        ++i;
        ++j;
      }
    }
  }
  else if (TA == 'C' && TB == 'C')
  {
    while (i < n1 && j < n2)
    {
      if (*(indx1 + i) < *(indx2 + j))
        ++i;
      else if (*(indx2 + j) < *(indx1 + i))
        ++j;
      else
      {
        res += ma::conj(*(A1 + i)) * ma::conj(*(A2 + j));
        ++i;
        ++j;
      }
    }
  }
  return res;
}

/*
template<class CSR,
         class MultiArray2D,
         typename = typename std::enable_if_t<(std::decay<CSR>::type::dimensionality == -2)>,
         typename = typename std::enable_if_t<(MultiArray2D::dimensionality==2)>
        >
void Matrix2MA(char TA, CSR const& A, MultiArray2D& M)
{
  using Type = typename MultiArray2D::element;
  using int_type = typename CSR::int_type;
  assert(TA=='N' || TA=='H' || TA=='T' || TA=='Z');
  if(TA=='N' || TA=='Z')
    M.reextent({A.size(0),A.size(1)});
  else if(TA=='T' || TA=='H')
    M.reextent({A.size(1),A.size(0)});
  else
    throw std::runtime_error(" Error: Unknown operation in Matrix2MA.\n");
  using std::fill_n;
  fill_n(M.origin(),M.num_elements(),Type(0));
  auto pbegin = A.pointers_begin();
  auto pend = A.pointers_end();
  int_type p0(pbegin[0]);
  auto v0 = A.non_zero_values_data();
  auto c0 = A.non_zero_indices2_data();
  if(TA=='N') {
    for(int i=0; i<A.size(0); i++)
      for(int_type ip=pbegin[i], ipend=pend[i]; ip<ipend; ip++)
        M[i][c0[ip-p0]] = Type(v0[ip-p0]);
  } else if(TA=='Z') {
    for(int i=0; i<A.size(0); i++)
      for(int_type ip=pbegin[i], ipend=pend[i]; ip<ipend; ip++)
        M[i][c0[ip-p0]] = ma::conj(Type(v0[ip-p0]));
  } else if(TA=='T') {
    for(int i=0; i<A.size(0); i++)
      for(int_type ip=pbegin[i], ipend=pend[i]; ip<ipend; ip++)
        M[c0[ip-p0]][i] = Type(v0[ip-p0]);
  } else if(TA=='H') {
    for(int i=0; i<A.size(0); i++)
      for(int_type ip=pbegin[i], ipend=pend[i]; ip<ipend; ip++)
        M[c0[ip-p0]][i] = ma::conj(Type(v0[ip-p0]));
  }
}

template<class CSR,
         class MultiArray2D,
         typename = typename std::enable_if_t<(std::decay<CSR>::type::dimensionality == -2)>,
         typename = typename std::enable_if_t<(MultiArray2D::dimensionality==2)>
        >
void Matrix2MAREF(char TA, CSR const& A, MultiArray2D& M)
{
  using Type = typename MultiArray2D::element;
  using int_type = typename CSR::int_type;
  assert(TA=='N' || TA=='H' || TA=='T' || TA=='Z');
  if( (TA=='N' || TA=='Z') && ( (M.size(0)!=A.size(0)) || (M.size(1)!=A.size(1)) ) )
    throw std::runtime_error(" Error: Wrong dimensions in Matrix2MAREF.\n");
  else if( (TA=='T' || TA=='H') && ( (M.size(0)!=A.size(1)) || (M.size(1)!=A.size(0)) ) )
    throw std::runtime_error(" Error: Wrong dimensions in Matrix2MAREF.\n");
  using std::fill_n;
  fill_n(M.origin(),M.num_elements(),Type(0));
  auto pbegin = A.pointers_begin();
  auto pend = A.pointers_end();
  int_type p0(pbegin[0]);
  auto v0 = A.non_zero_values_data();
  auto c0 = A.non_zero_indices2_data();
  if(TA=='N') {
    for(int i=0; i<A.size(0); i++)
      for(int_type ip=pbegin[i], ipend=pend[i]; ip<ipend; ip++)
        M[i][c0[ip-p0]] = Type(v0[ip-p0]);
  } else if(TA=='Z') {
    for(int i=0; i<A.size(0); i++)
      for(int_type ip=pbegin[i], ipend=pend[i]; ip<ipend; ip++)
        M[i][c0[ip-p0]] = ma::conj(Type(v0[ip-p0]));
  } else if(TA=='T') {
    for(int i=0; i<A.size(0); i++)
      for(int_type ip=pbegin[i], ipend=pend[i]; ip<ipend; ip++)
        M[c0[ip-p0]][i] = Type(v0[ip-p0]);
  } else if(TA=='H') {
    for(int i=0; i<A.size(0); i++)
      for(int_type ip=pbegin[i], ipend=pend[i]; ip<ipend; ip++)
        M[c0[ip-p0]][i] = ma::conj(Type(v0[ip-p0]));
  }
}
*/

/* Chooses rows of A based on occups vector and performs CSF2MA on subset of rows */
/*
template<class CSR,
         class MultiArray2D,
         class Vector,
         typename = typename std::enable_if_t<(std::decay<CSR>::type::dimensionality == -2)>,
         typename = typename std::enable_if_t<(MultiArray2D::dimensionality==2)>
        >
void Matrix2MA(char TA, CSR const& A, MultiArray2D& M, Vector const& occups)
{
  using Type = typename MultiArray2D::element;
  if(occups.size()==0) throw std::runtime_error(" Error: Empty occupation array in Matrix2MA.\n");
  assert(occups.size() <= A.size(0));
  int nrows = occups.size();
  assert(TA=='N' || TA=='H' || TA=='T' || TA=='Z');
  if(TA=='N' || TA=='Z') {
    if(M.size(0) != nrows || M.size(1) != A.size(1))
      M.reextent({nrows,A.size(1)});
  } else if(TA=='T' || TA=='H') {
    if(M.size(1) != nrows || M.size(0) != A.size(1))
      M.reextent({A.size(1),nrows});
  } else
    throw std::runtime_error(" Error: Unknown operation in Matrix2MA.\n");
  std::fill_n(M.origin(),M.num_elements(),Type(0));
  auto pbegin = A.pointers_begin();
  auto pend = A.pointers_end();
  auto p0 = pbegin[0];
  auto v0 = A.non_zero_values_data();
  auto c0 = A.non_zero_indices2_data();
  if(TA=='N') {
    for(int i=0; i<nrows; i++) {
      assert(occups[i] >= 0 && occups[i] < A.size(0));  
      int ik = occups[i];  
      for(int ip=pbegin[ik]; ip<pend[ik]; ip++)
        M[i][c0[ip-p0]] = static_cast<Type>(v0[ip-p0]);
    }    
  } else if(TA=='Z') {
    for(int i=0; i<nrows; i++) {
      assert(occups[i] >= 0 && occups[i] < A.size(0));
      int ik = occups[i];
      for(int ip=pbegin[ik]; ip<pend[ik]; ip++)
        M[i][c0[ip-p0]] = static_cast<Type>(ma::conj(v0[ip-p0]));
    }
  } else if(TA=='T') {
    for(int i=0; i<nrows; i++) {
      assert(occups[i] >= 0 && occups[i] < A.size(0));  
      int ik = occups[i];  
      for(int ip=pbegin[ik]; ip<pend[ik]; ip++)
        M[c0[ip-p0]][i] = static_cast<Type>(v0[ip-p0]);
    }    
  } else if(TA=='H') {
    for(int i=0; i<nrows; i++) {
      assert(occups[i] >= 0 && occups[i] < A.size(0));  
      int ik = occups[i];  
      for(int ip=pbegin[ik]; ip<pend[ik]; ip++)
        M[c0[ip-p0]][i] = static_cast<Type>(ma::conj(v0[ip-p0]));
    }    
  }
}
*/
namespace shm
{
template<class csr_matrix_out,
         class csr_matrix_in,
         typename = typename std::enable_if<std::decay<csr_matrix_in>::type::sparse>::type,
         typename = typename std::enable_if<std::decay<csr_matrix_in>::type::dimensionality == -2>::type,
         typename = typename std::enable_if<std::decay<csr_matrix_out>::type::sparse>::type,
         typename = typename std::enable_if<std::decay<csr_matrix_out>::type::dimensionality == -2>::type>
csr_matrix_out transpose(csr_matrix_in&& A)
{
  using integer    = typename std::decay<csr_matrix_in>::type::index_type;
  using value_type = typename std::decay<csr_matrix_out>::type::value_type;
  auto& comm       = *A.getAlloc().commP_;
  std::vector<std::size_t> sz_per_row(A.size(1));
  integer r0, rN, ncols = integer(A.size(1));
  integer rank = comm.rank(), size = comm.size();
  std::tie(r0, rN) = FairDivideBoundary(rank, ncols, size);

  auto pb0 = *(A.pointers_begin(0));
  for (integer i = 0; i < integer(A.size(0)); i++)
  {
    auto pbi = *(A.pointers_begin(i));
    auto pei = *(A.pointers_end(i));
    auto c0  = std::lower_bound(to_address((A.non_zero_indices2_data() + (pbi - pb0))),
                               to_address((A.non_zero_indices2_data() + (pei - pb0))), r0);
    auto cN  = std::lower_bound(to_address((A.non_zero_indices2_data() + (pbi - pb0))),
                               to_address((A.non_zero_indices2_data() + (pei - pb0))), rN);
    for (; c0 != cN; ++c0)
      ++sz_per_row[*c0];
  }
  comm.all_reduce_in_place_n(sz_per_row.begin(), sz_per_row.size(), std::plus<>());
  csr_matrix_out csr(std::tuple<std::size_t, std::size_t>{A.size(1), A.size(0)},
                     std::tuple<std::size_t, std::size_t>{0, 0}, sz_per_row, A.getAlloc());
  for (integer i = 0; i < integer(A.size(0)); i++)
  {
    auto pbi = *(A.pointers_begin(i));
    auto pei = *(A.pointers_end(i));
    auto c0  = std::lower_bound(to_address((A.non_zero_indices2_data() + (pbi - pb0))),
                               to_address((A.non_zero_indices2_data() + (pei - pb0))), r0);
    auto cN  = std::lower_bound(to_address((A.non_zero_indices2_data() + (pbi - pb0))),
                               to_address((A.non_zero_indices2_data() + (pei - pb0))), rN);
    auto dn  = std::distance(to_address(A.non_zero_indices2_data()), c0);
    auto v0  = A.non_zero_values_data() + dn;
    for (; c0 != cN; ++c0, ++v0)
      csr.emplace_back({*c0, i}, static_cast<value_type>(*v0));
  }
  comm.barrier();
  return csr;
}

template<class csr_matrix,
         class MultiArray2D,
         typename = typename std::enable_if<std::decay<csr_matrix>::type::sparse>::type,
         typename = typename std::enable_if<std::decay<csr_matrix>::type::dimensionality == -2>::type,
         typename = typename std::enable_if<std::decay<MultiArray2D>::type::dimensionality == 2>::type>
MultiArray2D transpose(csr_matrix&& A, MultiArray2D&& AT)
{
  using integer = typename std::decay<csr_matrix>::type::index_type;
  using Type    = typename std::decay<MultiArray2D>::type::element;
  assert(std::get<0>(A.sizes()) == std::get<1>(AT.sizes()));
  assert(std::get<1>(A.sizes()) == std::get<0>(AT.sizes()));
  auto& comm = *A.getAlloc().commP_;
  integer r0, rN, nrows = integer(A.size(0));
  integer rank = comm.rank(), size = comm.size();
  std::tie(r0, rN) = FairDivideBoundary(rank, nrows, size);

  auto pb0 = *(A.pointers_begin(0));
  for (integer i = r0; i < rN; i++)
  {
    auto pbi = *(A.pointers_begin(i));
    auto pei = *(A.pointers_end(i));
    auto c0  = A.non_zero_indices2_data() + (pbi - pb0);
    auto v0  = A.non_zero_values_data() + (pbi - pb0);
    for (; pbi != pei; ++pbi, ++c0, ++v0)
      AT[*c0][i] = static_cast<Type>(*v0);
  }
  comm.barrier();
  return std::forward<MultiArray2D>(AT);
}

} // namespace shm
} // namespace csr
#endif
