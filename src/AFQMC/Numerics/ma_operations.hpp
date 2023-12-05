//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Alfredo Correa, correaa@llnl.gov
//    Lawrence Livermore National Laboratory
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Alfredo Correa, correaa@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef MA_OPERATIONS_HPP
#define MA_OPERATIONS_HPP

#include "ma_blas.hpp"
#include "ma_lapack.hpp"
#include "AFQMC/Numerics/detail/sparse.hpp"
#include "AFQMC/Numerics/determinant.hpp"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"

#include <utility>     // declval
#include <type_traits> // enable_if
#include <vector>

#include "AFQMC/Utilities/type_conversion.hpp"

// pass all pointers through a dispatch() call that can be customized to map pointer types
// then define it generically to returned the argument, and specialize it for array_ptr_raw_ptr_dispatch<T>

namespace ma
{
using qmcplusplus::afqmc::pointedType;
using qmcplusplus::afqmc::to_address;

template<class MultiArray2D, typename = typename std::enable_if<(MultiArray2D::dimensionality > 1)>::type>
bool is_hermitian(MultiArray2D const& A)
{
  using ma::conj;
  if (A.size() != std::get<1>(A.sizes()))
    return false;
  for (int i = 0; i != std::get<0>(A.sizes()); ++i)
    for (int j = i + 1; j != std::get<1>(A.sizes()); ++j)
      if (std::abs(A[i][j] - ma::conj(A[j][i])) > 1e-12)
        return false;
  return true;
}

template<class MultiArray2D, typename = typename std::enable_if<(MultiArray2D::dimensionality > 1)>::type>
bool is_symmetric(MultiArray2D const& A)
{
  if (A.size() != A.size(1))
    return false;
  for (int i = 0; i != A.size(); ++i)
    for (int j = i + 1; j != A.size(1); ++j)
      if (std::abs(A[i][j] - A[j][i]) > 1e-12)
        return false;
  return true;
}

template<class MA>
struct op_tag : std::integral_constant<char, 'N'>
{}; // see specializations
template<class MA>
MA arg(MA&& ma)
{
  return std::forward<MA>(ma);
} // see specializations below

template<class MultiArray2D,
         typename = typename std::enable_if<(std::decay<MultiArray2D>::type::dimensionality > 1)>::type>
MultiArray2D&& transpose(MultiArray2D&& A)
{
  assert(A.size() == A.size(1));
  typename std::decay<MultiArray2D>::type::element one(1.0);
  return ma::geam<'T'>(one, arg(A), std::forward<MultiArray2D>(A));
}

template<class MultiArray2D,
         typename = typename std::enable_if<(std::decay<MultiArray2D>::type::dimensionality == 2)>::type>
void transform(MultiArray2D&& A)
{
  assert(arg(A).size() == arg(A).size(1));
  typename std::decay<MultiArray2D>::type::element one(1.0);
  ma::geam<op_tag<typename std::decay<MultiArray2D>::type>::value>(one, arg(A), arg(A));
}

template<class MultiArray2DA,
         class MultiArray2DB,
         typename = typename std::enable_if<(std::decay<MultiArray2DA>::type::dimensionality > 1)>::type,
         typename = typename std::enable_if<(std::decay<MultiArray2DB>::type::dimensionality > 1)>::type,
         typename = typename std::enable_if<(std::decay<MultiArray2DA>::type::dimensionality ==
                                             std::decay<MultiArray2DB>::type::dimensionality)>::type>
MultiArray2DB&& transpose(MultiArray2DA&& A, MultiArray2DB&& B)
{
  typename std::decay<MultiArray2DA>::type::element one(1.0);
  return ma::geam<'T'>(one, arg(A), std::forward<MultiArray2DB>(B));
}

template<class MultiArray2DA,
         class MultiArray2DB,
         typename = typename std::enable_if<(std::decay<MultiArray2DA>::type::dimensionality > 1)>::type,
         typename = typename std::enable_if<(std::decay<MultiArray2DB>::type::dimensionality > 1)>::type,
         typename = typename std::enable_if<(std::decay<MultiArray2DA>::type::dimensionality ==
                                             std::decay<MultiArray2DB>::type::dimensionality)>::type>
MultiArray2DB&& transform(MultiArray2DA const& A, MultiArray2DB&& B)
{
  typename std::decay<MultiArray2DA>::type::element one(1.0);
  return ma::geam<op_tag<MultiArray2DA>::value>(one, arg(A), std::forward<MultiArray2DB>(B));
}

template<
    class T,
    class MultiArray2DA,
    class MultiArray2DB,
    class MultiArray2DC,
    typename = typename std::enable_if<MultiArray2DA::dimensionality == 2 and MultiArray2DB::dimensionality == 2 and
                                       std::decay<MultiArray2DC>::type::dimensionality == 2>::type>
MultiArray2DC&& add(T alpha, MultiArray2DA const& A, T beta, MultiArray2DB const& B, MultiArray2DC&& C)
{
  return ma::geam<op_tag<MultiArray2DA>::value, op_tag<MultiArray2DB>::value>(alpha, arg(A), beta, arg(B),
                                                                              std::forward<MultiArray2DC>(C));
}

template<
    class T,
    class MultiArray3DA,
    class MultiArray3DB,
    class MultiArray3DC,
    typename = typename std::enable_if<MultiArray3DA::dimensionality == 3 and MultiArray3DB::dimensionality == 3 and
                                       std::decay<MultiArray3DC>::type::dimensionality == 3>::type>
MultiArray3DC&& productStridedBatched(T alpha,
                                      MultiArray3DA const& A,
                                      MultiArray3DB const& B,
                                      T beta,
                                      MultiArray3DC&& C)
{
  return ma::gemmStridedBatched<op_tag<MultiArray3DB>::value, op_tag<MultiArray3DA>::value>(alpha, arg(B), arg(A), beta,
                                                                                            std::forward<MultiArray3DC>(
                                                                                                C));
}

template<
    class T,
    class MultiArray2DA,
    class MultiArray1DB,
    class MultiArray1DC,
    typename = typename std::enable_if<MultiArray2DA::dimensionality == 2 and MultiArray1DB::dimensionality == 1 and
                                       std::decay<MultiArray1DC>::type::dimensionality == 1>::type,
    typename = void // TODO change to use dispatch
    >
MultiArray1DC&& product(T alpha, MultiArray2DA const& A, MultiArray1DB const& B, T beta, MultiArray1DC&& C)
{
  return ma::gemv < (op_tag<MultiArray2DA>::value == 'N')
      ? 'T'
      : 'N' > (alpha, arg(A), B, beta, std::forward<MultiArray1DC>(C));
}

// sparse matrix-MultiArray<1> interface
template<
    class T,
    class SparseMatrixA,
    class MultiArray1DB,
    class MultiArray1DC,
    typename = typename std::enable_if<SparseMatrixA::dimensionality == -2 and MultiArray1DB::dimensionality == 1 and
                                       std::decay<MultiArray1DC>::type::dimensionality == 1>::type,
    typename = void,
    typename = void, // TODO change to use dispatch
    typename = void  // TODO change to use dispatch
    >
MultiArray1DC&& product(T alpha, SparseMatrixA const& A, MultiArray1DB const& B, T beta, MultiArray1DC&& C)
{
  using elementA = typename SparseMatrixA::element;
  using elementB = typename MultiArray1DB::element;
  using elementC = typename std::decay<MultiArray1DC>::type::element;
  //static_assert(std::is_same<elementA,elementB>::value,"Problems with sparse dispatch");
  //static_assert(std::is_same<elementA,elementC>::value,"Problems with sparse dispatch");
  //static_assert(std::is_convertible<T,elementC>::value,"Problems with sparse dispatch");
  assert(op_tag<SparseMatrixA>::value == 'N' || op_tag<SparseMatrixA>::value == 'T');
  assert(op_tag<MultiArray1DB>::value == 'N');
  if (op_tag<SparseMatrixA>::value == 'N')
  {
    assert(arg(A).size() == std::forward<MultiArray1DC>(C).size());
    assert(arg(A).size(1) == arg(B).size());
  }
  else
  {
    assert(arg(A).size() == arg(B).size());
    assert(arg(A).size(1) == std::forward<MultiArray1DC>(C).size());
  }

  csrmv(op_tag<SparseMatrixA>::value, arg(A).size(), arg(A).size(1), elementA(alpha), "GxxCxx",
        pointer_dispatch(arg(A).non_zero_values_data()), pointer_dispatch(arg(A).non_zero_indices2_data()),
        pointer_dispatch(arg(A).pointers_begin()), pointer_dispatch(arg(A).pointers_end()),
        pointer_dispatch(arg(B).origin()), elementA(beta), pointer_dispatch(C.origin()));

  return std::forward<MultiArray1DC>(C);
}

///*
template<class MultiArray2DA,
         class MultiArray1DB,
         class MultiArray1DC,
         typename = typename std::enable_if<
             (MultiArray2DA::dimensionality == 2 or MultiArray2DA::dimensionality == -2) and
             MultiArray1DB::dimensionality == 1 and std::decay<MultiArray1DC>::type::dimensionality == 1>::type,
         typename = void // TODO change to use dispatch
         >
MultiArray1DC&& product(MultiArray2DA const& A, MultiArray1DB const& B, MultiArray1DC&& C)
{
  using Type = typename std::decay<MultiArray2DA>::type::element;
  return product(Type(1.), A, B, Type(0.), std::forward<MultiArray1DC>(C));
}
//*/

template<
    class T,
    class MultiArray2DA,
    class MultiArray2DB,
    class MultiArray2DC,
    typename = typename std::enable_if<MultiArray2DA::dimensionality == 2 and MultiArray2DB::dimensionality == 2 and
                                       std::decay<MultiArray2DC>::type::dimensionality == 2>::type>
MultiArray2DC&& product(T alpha, MultiArray2DA const& A, MultiArray2DB const& B, T beta, MultiArray2DC&& C)
{
  return ma::gemm<op_tag<MultiArray2DB>::value, op_tag<MultiArray2DA>::value>(alpha, arg(B), arg(A), beta,
                                                                              std::forward<MultiArray2DC>(C));
}

// sparse matrix-MultiArray interface
template<
    class T,
    class SparseMatrixA,
    class MultiArray2DB,
    class MultiArray2DC,
    typename = typename std::enable_if<SparseMatrixA::dimensionality == -2 and MultiArray2DB::dimensionality == 2 and
                                       std::decay<MultiArray2DC>::type::dimensionality == 2>::type,
    typename = void,
    typename = void // TODO change to use dispatch
    >
MultiArray2DC&& product(T alpha, SparseMatrixA const& A, MultiArray2DB const& B, T beta, MultiArray2DC&& C)
{
  using elementA = std::remove_cv_t<typename SparseMatrixA::element>;
  using elementB = std::remove_cv_t<typename MultiArray2DB::element>;
  using elementC = typename std::decay<MultiArray2DC>::type::element;
  //static_assert(std::is_same<elementA,elementB>::value,"Problems with sparse dispatch");
  //static_assert(std::is_same<elementA,elementC>::value,"Problems with sparse dispatch");
  //static_assert(std::is_convertible<T,elementC>::value,"Problems with sparse dispatch");
  assert(op_tag<SparseMatrixA>::value == 'N' || op_tag<SparseMatrixA>::value == 'T' ||
         op_tag<SparseMatrixA>::value == 'C' || op_tag<SparseMatrixA>::value == 'H');
  assert(op_tag<MultiArray2DB>::value == 'N');
  assert(arg(B).stride(1) == 1);
  assert(std::forward<MultiArray2DC>(C).stride(1) == 1);
  if (op_tag<SparseMatrixA>::value == 'N')
  {
    assert(arg(A).size() == std::forward<MultiArray2DC>(C).size());
    assert( arg(A).size(1) == arg(B).size() );
    assert( std::get<1>(arg(B).sizes()) == std::get<1>(std::forward<MultiArray2DC>(C).sizes()) );
  }
  else
  {
    assert(arg(A).size() == arg(B).size());
    assert(std::get<1>(arg(A).sizes()) == std::forward<MultiArray2DC>(C).size());
    assert(std::get<1>(arg(B).sizes()) == std::get<1>(std::forward<MultiArray2DC>(C).sizes()));
  }

  csrmm(op_tag<SparseMatrixA>::value, arg(A).size(), std::get<1>(arg(B).sizes()), std::get<1>(arg(A).sizes()), elementA(alpha), "GxxCxx",
        pointer_dispatch(arg(A).non_zero_values_data()), pointer_dispatch(arg(A).non_zero_indices2_data()),
        pointer_dispatch(arg(A).pointers_begin()), pointer_dispatch(arg(A).pointers_end()),
        pointer_dispatch(arg(B).origin()), arg(B).stride(), elementA(beta), pointer_dispatch(C.origin()), C.stride());

  return std::forward<MultiArray2DC>(C);
}

template<
    class T,
    class SparseMatrixA,
    class MultiArray2DB,
    class MultiArray2DC,
    typename = typename std::enable_if<SparseMatrixA::dimensionality == -2 and MultiArray2DB::dimensionality == 2 and
                                       std::decay<MultiArray2DC>::type::dimensionality == 2>::type,
    typename = void,
    typename = void // TODO change to use dispatch
    >
MultiArray2DC&& product(T alpha, MultiArray2DB const& B, SparseMatrixA const& A, T beta, MultiArray2DC&& C)
{
  APP_ABORT(" Error: DenseMatrix x SparseMatrix product not implemented. \n\n\n");
  return std::forward<MultiArray2DC>(C);
}

template<
    class MultiArray2DA,
    class MultiArray2DB,
    class MultiArray2DC,
    typename = typename std::enable_if<(MultiArray2DA::dimensionality == 2 or MultiArray2DA::dimensionality == -2) and
                                       (MultiArray2DB::dimensionality == 2 or MultiArray2DB::dimensionality == -2) and
                                       std::decay<MultiArray2DC>::type::dimensionality == 2>::type>
MultiArray2DC&& product(MultiArray2DA const& A, MultiArray2DB const& B, MultiArray2DC&& C)
{
  using Type = typename std::decay<MultiArray2DA>::type::element;
  return product(Type(1.), A, B, Type(0.), std::forward<MultiArray2DC>(C));
}

template<class T,
         class MultiArrayPtr2DA,
         class MultiArrayPtr2DB,
         class MultiArrayPtr2DC,
         typename = typename std::enable_if<(pointedType<MultiArrayPtr2DA>::dimensionality == 2) and
                                            (pointedType<MultiArrayPtr2DB>::dimensionality == 2) and
                                            pointedType<MultiArrayPtr2DC>::dimensionality == 2>::type,
         typename = void,
         typename = void>
void BatchedProduct(char TA,
                    char TB,
                    T alpha,
                    std::vector<MultiArrayPtr2DA>& A,
                    std::vector<MultiArrayPtr2DB>& B,
                    T beta,
                    std::vector<MultiArrayPtr2DC>& C)
{
  int nbatch = C.size();
  assert(A.size() >= nbatch);
  assert(B.size() >= nbatch);

  // need to do it this way to preserve qualifiers! can also use decltype
  //using pointerA_ = typename pointedType<MultiArrayPtr2DA>::element_ptr;
  //using pointerB_ = typename pointedType<MultiArrayPtr2DB>::element_ptr;
  //using pointerC_ = typename pointedType<MultiArrayPtr2DC>::element_ptr;
  using pointerA = decltype(pointer_dispatch((*A[0]).origin()));
  using pointerB = decltype(pointer_dispatch((*B[0]).origin()));
  using pointerC = decltype(pointer_dispatch((*C[0]).origin()));
  using element  = typename pointedType<MultiArrayPtr2DA>::element;

  int M = std::get<1>((*C[0]).sizes());
  int N = (*C[0]).size();
  int K;
  if (TB == 'N')
    K = (*B[0]).size();
  else
    K = std::get<1>((*B[0]).sizes());
  int lda = (*A[0]).stride();
  int ldb = (*B[0]).stride();
  int ldc = (*C[0]).stride();
  std::vector<pointerA> Ai;
  std::vector<pointerB> Bi;
  std::vector<pointerC> Ci;
  Ai.reserve(nbatch);
  Bi.reserve(nbatch);
  Ci.reserve(nbatch);
  for (int i = 0; i < nbatch; i++)
  {
    assert(lda == (*A[i]).stride());
    assert(ldb == (*B[i]).stride());
    assert(ldc == (*C[i]).stride());
    assert(M == std::get<1>((*C[i]).sizes()));
    assert(N == (*C[i]).size());
    if (TB == 'N')
    {
      assert(K == (*B[i]).size());
      assert(M == std::get<1>((*B[i]).sizes()));
    }
    else
    {
      assert(K == std::get<1>((*B[i]).sizes()));
      assert(M == (*B[i]).size());
    }
    if (TA == 'N')
    {
      assert(K == std::get<1>((*A[i]).sizes()));
      assert(N == (*A[i]).size());
    }
    else
    {
      assert(K == (*A[i]).size());
      assert(N == std::get<1>((*A[i]).sizes()));
    }
    Ai.emplace_back(pointer_dispatch((*A[i]).origin()));
    Bi.emplace_back(pointer_dispatch((*B[i]).origin()));
    Ci.emplace_back(pointer_dispatch((*C[i]).origin()));
  }

  using ma::gemmBatched;
  gemmBatched(TB, TA, M, N, K, element(alpha), Bi.data(), ldb, Ai.data(), lda, element(beta), Ci.data(), ldc, nbatch);
}

// no batched sparse product yet, serialize call
template<class T,
         class MultiArrayPtr2DA,
         class MultiArrayPtr2DB,
         class MultiArrayPtr2DC,
         typename = typename std::enable_if<(pointedType<MultiArrayPtr2DA>::dimensionality == -2) and
                                            (pointedType<MultiArrayPtr2DB>::dimensionality == 2) and
                                            pointedType<MultiArrayPtr2DC>::dimensionality == 2>::type,
         typename = void>
void BatchedProduct(char TA,
                    char TB,
                    T alpha,
                    std::vector<MultiArrayPtr2DA>& A,
                    std::vector<MultiArrayPtr2DB>& B,
                    T beta,
                    std::vector<MultiArrayPtr2DC>& C)
{
  int nbatch = C.size();
  assert(A.size() >= nbatch);
  assert(B.size() >= nbatch);

  using elementA = std::remove_cv_t<typename pointedType<MultiArrayPtr2DA>::element>;
  using elementB = std::remove_cv_t<typename pointedType<MultiArrayPtr2DB>::element>;
  using elementC = typename pointedType<MultiArrayPtr2DC>::element;
  assert(TA == 'N' || TA == 'T' || TA == 'C' || TA == 'H');
  assert(TB == 'N');
  /*
        if(op_tag<SparseMatrixA>::value == 'N') {
            assert(arg(A).size(0) == std::forward<MultiArray2DC>(C).size(0));
            assert(arg(A).size(1) == arg(B).size(0));
            assert(arg(B).size(1) == std::forward<MultiArray2DC>(C).size(1));
        } else {
            assert(arg(A).size(0) == arg(B).size(0));
            assert(arg(A).size(1) == std::forward<MultiArray2DC>(C).size(0));
            assert(arg(B).size(1) == std::forward<MultiArray2DC>(C).size(1));
        }
*/

  for (int i = 0; i < nbatch; i++)
  {
    csrmm(TA, (*A[i]).size(), std::get<1>((*B[i]).sizes()), std::get<1>((*A[i]).sizes()), elementA(alpha), "GxxCxx",
          pointer_dispatch((*A[i]).non_zero_values_data()), pointer_dispatch((*A[i]).non_zero_indices2_data()),
          pointer_dispatch((*A[i]).pointers_begin()), pointer_dispatch((*A[i]).pointers_end()),
          pointer_dispatch((*B[i]).origin()), (*B[i]).stride(), elementA(beta), pointer_dispatch((*C[i]).origin()),
          (*C[i]).stride());
  }
}

template<class MultiArrayPtr2DA,
         class MultiArrayPtr2DB,
         class MultiArrayPtr2DC,
         typename = typename std::enable_if<(pointedType<MultiArrayPtr2DA>::dimensionality == 2 or
                                             pointedType<MultiArrayPtr2DA>::dimensionality == -2) and
                                            (pointedType<MultiArrayPtr2DB>::dimensionality == 2) and
                                            pointedType<MultiArrayPtr2DC>::dimensionality == 2>::type>
void BatchedProduct(char TA,
                    char TB,
                    std::vector<MultiArrayPtr2DA>& A,
                    std::vector<MultiArrayPtr2DB>& B,
                    std::vector<MultiArrayPtr2DC>& C)
{
  using Type = typename pointedType<MultiArrayPtr2DA>::element;
  return BatchedProduct(TA, TB, Type(1.), A, B, Type(0.), C);
}


template<class MultiArray2D>
struct normal_tag
{
  MultiArray2D arg1;
  static auto const dimensionality = std::decay<MultiArray2D>::type::dimensionality;
  using element                    = typename std::decay<MultiArray2D>::type::element;
  normal_tag(normal_tag const&)    = delete;
  normal_tag(normal_tag&&)         = default;
  static const char tag            = 'N';
};

template<class MultiArray2D>
struct op_tag<normal_tag<MultiArray2D>> : std::integral_constant<char, 'N'>
{};

template<class MultiArray2D>
normal_tag<MultiArray2D> normal(MultiArray2D&& arg)
{
  return {std::forward<MultiArray2D>(arg)};
}

template<class MultiArray2D>
MultiArray2D&& arg(normal_tag<MultiArray2D> const& nt)
{
  return nt.arg1;
}

template<class MultiArray2D>
struct transpose_tag
{
  MultiArray2D arg1;
  static auto const dimensionality    = std::decay<MultiArray2D>::type::dimensionality;
  using element                       = typename std::decay<MultiArray2D>::type::element;
  transpose_tag(transpose_tag const&) = delete;
  transpose_tag(transpose_tag&&)      = default;
  static const char tag               = 'T';
};

template<class MultiArray2D>
struct op_tag<transpose_tag<MultiArray2D>> : std::integral_constant<char, 'T'>
{};

template<class MultiArray2D>
transpose_tag<MultiArray2D&&> transposed_matrix(MultiArray2D&& arg)
{
  return {std::forward<MultiArray2D>(arg)};
}

// return a pointer instead of a copy, that way you don't care if it is copyable or not
template<class MultiArray2D>
MultiArray2D const& arg(transpose_tag<MultiArray2D> const& tt)
{
  return tt.arg1;
}

template<class MultiArray2D>
struct hermitian_tag
{
  MultiArray2D arg1;
  static auto const dimensionality    = std::decay<MultiArray2D>::type::dimensionality;
  using element                       = typename std::decay<MultiArray2D>::type::element;
  hermitian_tag(hermitian_tag const&) = delete;
  hermitian_tag(hermitian_tag&&)      = default;
  static const char tag               = 'C';
};

template<class MultiArray2D>
hermitian_tag<MultiArray2D&&> hermitian(MultiArray2D&& arg)
{
  return {std::forward<MultiArray2D>(arg)};
}

template<class MultiArray2D>
MultiArray2D const& arg(hermitian_tag<MultiArray2D> const& nt)
{
  return nt.arg1;
}

template<class MultiArray2D>
struct op_tag<hermitian_tag<MultiArray2D>> : std::integral_constant<char, 'C'>
{};

//template<class MultiArray2D>
//MultiArray2D const& arg(hermitian_tag<MultiArray2D>& ht){return ht.arg1;}


template<class MA2D>
auto T(MA2D&& arg) -> decltype(transposed_matrix(std::forward<MA2D>(arg)))
{
  return transposed_matrix(std::forward<MA2D>(arg));
}
template<class MA2D>
auto H(MA2D&& arg) -> decltype(hermitian(std::forward<MA2D>(arg)))
{
  return hermitian(std::forward<MA2D>(arg));
}
template<class MA2D>
auto N(MA2D&& arg) -> decltype(normal(std::forward<MA2D>(arg)))
{
  return normal(std::forward<MA2D>(arg));
}

template<class MA2D>
auto trans(MA2D&& arg) -> decltype(transposed_matrix(std::forward<MA2D>(arg)))
{
  return transposed_matrix(std::forward<MA2D>(arg));
}
template<class MA2D>
auto herm(MA2D&& arg) -> decltype(hermitian(std::forward<MA2D>(arg)))
{
  return hermitian(std::forward<MA2D>(arg));
}
//template<class MA2D> auto norm(MA2D&& arg)
//->decltype(normal(std::forward<MA2D>(arg))){
//	return normal(std::forward<MA2D>(arg));
//}

template<class MultiArray2D>
int invert_optimal_workspace_size(MultiArray2D&& m)
{
  return std::max(getri_optimal_workspace_size(m), getrf_optimal_workspace_size(m));
}

template<class T, class MultiArray2D>
T invert(MultiArray2D&& m, T LogOverlapFactor)
{
  using element         = typename std::decay<MultiArray2D>::type::element;
  using allocator_type  = typename std::decay<MultiArray2D>::type::allocator_type;
  using iallocator_type = typename allocator_type::template rebind<int>::other;
  using extensions      = typename boost::multi::layout_t<1u>::extensions_type;
  using qmcplusplus::afqmc::fill2D;
  auto bufferSize(invert_optimal_workspace_size(std::forward<MultiArray2D>(m)));
  boost::multi::array<element, 1, allocator_type> WORK(extensions{bufferSize}, m.get_allocator());
  boost::multi::array<int, 1, iallocator_type> pivot(extensions{m.size() + 1}, iallocator_type{m.get_allocator()});

  getrf(std::forward<MultiArray2D>(m), pivot, WORK);
  T detvalue = determinant_from_getrf<T>(m.size(), pointer_dispatch(m.origin()), m.stride(),
                                         pointer_dispatch(pivot.data()), LogOverlapFactor);
  if (std::abs(detvalue) == 0.0)
    fill2D(m.size(), std::get<1>(m.sizes()), pointer_dispatch(m.origin()), m.stride(), element(0.0));
  else
    getri(std::forward<MultiArray2D>(m), pivot, WORK);
  return detvalue;
}

template<class T, class MultiArray2D, class MultiArray1D, class Buffer>
T invert(MultiArray2D&& m, MultiArray1D&& pivot, Buffer&& WORK, T LogOverlapFactor)
{
  assert(m.size() == std::get<1>(m.sizes()));
  assert(pivot.size() >= m.size() + 1);
  using element = typename std::decay<MultiArray2D>::type::element;
  using qmcplusplus::afqmc::fill2D;

  getrf(std::forward<MultiArray2D>(m), pivot, WORK);
  T detvalue = determinant_from_getrf<T>(m.size(), pointer_dispatch(m.origin()), m.stride(),
                                         pointer_dispatch(pivot.data()), LogOverlapFactor);
  if (std::abs(detvalue) == 0.0)
    fill2D(m.size(), std::get<1>(m.sizes()), pointer_dispatch(m.origin()), m.stride(), element(0.0));
  else
    getri(std::forward<MultiArray2D>(m), pivot, WORK);
  return detvalue;
}

/*
template<class MultiArray2D, 
         class MultiArray1DS, 
         class MultiArray2DU, 
         class MultiArray2DVT, 
         class Buffer, 
         class T = typename std::decay<MultiArray2D>::type::element>
void invert_withSVD(MultiArray2D&& m, MultiArray1DS&& S, MultiArray2DU&& U, MultiArray2DVT&& VT, Buffer&& WORK, T* detvalue){
        assert(m.size(0) == m.size(1));
        assert(S.size(0) == m.size(0));
        assert(VT.size(0) == m.size(0));
        assert(VT.size(0) == VT.size(1));
        using element = typename std::decay<MultiArray2D>::type::element;

        // m = U * S * VT
        // inv(m) = H(VT) * inv(S) * H(U)
        gesvd('A','A',std::forward<MultiArray2D>(m),S,U,VT,WORK);
        gesvd_determinant_and_regularization_of_singular_values<T>(m.size(0), 
                    pointer_dispatch(S.origin()), detvalue);
        // VT = VT * inv(S), which works since S is diagonal and real
        term_by_term_matrix_vector(TOp_DIV,1,VT.size(0),VT.size(1),pointer_dispatch(VT.origin()),
                    VT.stride(0),pointer_dispatch(S.origin()),1);
        product(H(VT),H(U),std::forward<MultiArray2D>(m));
}
*/

template<class T, class MultiArray2D, class MultiArray1D, class Buffer>
T determinant(MultiArray2D&& m, MultiArray1D&& pivot, Buffer&& WORK, T LogOverlapFactor)
{
  assert(m.size() == std::get<1>(m.sizes()));
  assert(pivot.size() >= m.size());

  getrf(std::forward<MultiArray2D>(m), std::forward<MultiArray1D>(pivot), WORK);
  return determinant_from_getrf<T>(m.size(), pointer_dispatch(m.origin()), m.stride(), pointer_dispatch(pivot.data()),
                                   LogOverlapFactor);
}

template<class MultiArray2D, typename = typename std::enable_if_t<MultiArray2D::dimensionality == 2>>
MultiArray2D exp(MultiArray2D const& A, bool printeV = false)
{
  using Type     = typename MultiArray2D::element;
  using RealType = typename qmcplusplus::afqmc::remove_complex<Type>::value_type;
  using TVec     = boost::multi::array<RealType, 1>;
  using TMat     = boost::multi::array<Type, 2>;
  using eigSys   = std::pair<TVec, TMat>;
  assert(A.size() == std::get<1>(A.sizes()));
  typename MultiArray2D::size_type N = A.size();

  MultiArray2D ExpA({N, N});
  std::fill_n(pointer_dispatch(ExpA.origin()), N * N, Type(0));

  if (is_hermitian(A))
  {
    // A = V*M*H(Z)
    eigSys V = symEig<TVec, TMat>(A);

    // exp(A) = V*exp(M)*H(Z)
    MultiArray2D TA({N, N});
    if (printeV)
    {
      qmcplusplus::app_log()
          << "***********************  Eigenvalues of exponentiated matrix *********************** \n";
      qmcplusplus::app_log() << " i    eigV[i]    exp(eigV[i]) \n";
      for (int j = 0; j < N; j++)
        qmcplusplus::app_log() << std::setprecision(16) << j << " " << V.first[j] << " " << std::exp(V.first[j])
                               << "\n";
      qmcplusplus::app_log()
          << "************************************************************************************ \n";
      qmcplusplus::app_log() << "\n\n" << std::endl;
    }
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        TA[i][j] = V.second[i][j] * std::exp(V.first[j]);
    product(TA, H(V.second), ExpA);
  }
  else
  {
    throw std::runtime_error("exp(A) not implemented for non-symmetric A");
  }

  return ExpA;
}

template<class MultiArray2D, typename = typename std::enable_if_t<MultiArray2D::dimensionality == 2>>
MultiArray2D&& cholesky(MultiArray2D&& A)
{
  assert(A.size() == std::get<1>(A.sizes()));
  if (is_hermitian(A))
  {
    return potrf(transpose(std::forward<MultiArray2D>(A)));
  }
  else
  {
    throw std::runtime_error("cholesky(A) not implemented for non-symmetric A");
    return std::forward<MultiArray2D>(A);
  }
}

template<class MultiArray2DA, class MultiArray2DB, class T>
bool equal(MultiArray2DB const& a, MultiArray2DA const& b, T tol = 0)
{
  if (a.size() != b.size() or a.size() != b.size())
    return false;
  using std::abs;
  for (int i = 0; i != a.size(); ++i)
    for (int j = 0; j != std::get<1>(a.sizes()); ++j)
      if (abs(a[i][j] - b[i][j]) > tol)
        return false;
  return true;
}

} // namespace ma


#ifdef _TEST_MA_OPERATIONS


#include <vector>
#include <iostream>

using std::cout;

int main()
{
  {
    std::vector<double> m = {
        9.,  24., 30., 4., 10.,
        12., 14., 16., 36. //,
                           //	9., 6., 1.
    };
    boost::multi::array_ref<double, 2> M(m.data(), {3, 3});
    assert(M.num_elements() == m.size());
    std::vector<double> x = {1., 2., 3.};
    boost::multi::array_ref<double, 1> X(x.data(), boost::iextensions<1u>{x.size()});
    std::vector<double> y(3);
    boost::multi::array_ref<double, 1> Y(y.data(), boost::iextensions<1u>{y.size()});

    using ma::T;
    ma::product(M, X, Y); // Y := M X

    std::vector<double> mx = {147., 60., 154.};
    boost::multi::array_ref<double, 1> MX(mx.data(), boost::iextensions<1u>{mx.size()});
    assert(MX == Y);
  }
  {
    std::vector<double> m = {9., 24., 30., 2., 4., 10., 12., 1., 14., 16., 36., 20.};
    boost::multi::array_ref<double, 2> M(m.data(), {3, 4});
    assert(M.num_elements() == m.size());
    std::vector<double> x = {1., 2., 3., 4.};
    boost::multi::array_ref<double, 1> X(x.data(), boost::iextensions<1u>{x.size()});
    std::vector<double> y(3);
    boost::multi::array_ref<double, 1> Y(y.data(), boost::iextensions<1u>{y.size()});

    using ma::T;
    ma::product(M, X, Y); // Y := M X

    std::vector<double> mx = {155., 64., 234.};
    boost::multi::array_ref<double, 1> MX(mx.data(), boost::iextensions<1u>{mx.size()});
    assert(MX == Y);
  }
  {
    std::vector<double> m = {9., 24., 30., 2., 4., 10., 12., 1., 14., 16., 36., 20.};
    boost::multi::array_ref<double, 2> M(m.data(), {3, 4});
    assert(M.num_elements() == m.size());
    std::vector<double> x = {1., 2., 3.};
    boost::multi::array_ref<double, 1> X(x.data(), boost::iextensions<1u>{x.size()});
    std::vector<double> y(4);
    boost::multi::array_ref<double, 1> Y(y.data(), boost::iextensions<1u>{y.size()});

    using ma::T;
    ma::product(T(M), X, Y); // Y := T(M) X

    std::vector<double> mx = {59., 92., 162., 64.};
    boost::multi::array_ref<double, 1> MX(mx.data(), boost::iextensions<1u>{mx.size()});
    assert(MX == Y);
  }
  {
    std::vector<double> m = {9., 24., 30., 9., 4., 10., 12., 7., 14., 16., 36., 1.};
    boost::multi::array_ref<double, 2> M(m.data(), {3, 4});
    std::vector<double> x = {1., 2., 3., 4.};
    boost::multi::array_ref<double, 1> X(x.data(), boost::iextensions<1u>{x.size()});
    std::vector<double> y = {4., 5., 6.};
    boost::multi::array_ref<double, 1> Y(y.data(), boost::iextensions<1u>{y.size()});
    ma::product(M, X, Y); // y := M x

    std::vector<double> y2 = {183., 88., 158.};
    boost::multi::array_ref<double, 1> Y2(y2.data(), boost::iextensions<1u>{y2.size()});
    assert(Y == Y2);
  }

  {
    std::vector<double> m = {1., 2., 1., 2., 5., 8., 1., 8., 9.};
    boost::multi::array_ref<double, 2> M(m.data(), {3, 3});
    assert(ma::is_hermitian(M));
  }
  {
    std::vector<double> m = {
        1., 0., 2., 0., 1., 0., 2., 0., 5., 0., 8., -1., 1., 0., 8., 1., 9., 0.,
    };
    boost::multi::array_ref<std::complex<double>, 2> M(reinterpret_cast<std::complex<double>*>(m.data()), {3, 3});
    assert(ma::is_hermitian(M));
  }
  {
    std::vector<double> m = {1., 2., 1., 2., 5., 8., 1., 8., 9.};
    boost::multi::array_ref<double, 2> M(m.data(), {3, 3});
    assert(ma::is_hermitian(M));
  }
  {
    std::vector<double> a = {1., 0., 1., 3., 5., 8., 4., 8., 9.};
    boost::multi::array_ref<double, 2> A(a.data(), {3, 3});
    assert(A.num_elements() == a.size());
    std::vector<double> b = {6., 2., 8., 9., 5., 5., 1., 7., 9.};
    boost::multi::array_ref<double, 2> B(b.data(), {3, 3});
    assert(B.num_elements() == b.size());

    std::vector<double> c(9);
    boost::multi::array_ref<double, 2> C(c.data(), {3, 3});
    assert(C.num_elements() == c.size());

    ma::product(A, B, C);

    std::vector<double> ab = {7., 9., 17., 71., 87., 121., 105., 111., 153.};
    boost::multi::array_ref<double, 2> AB(ab.data(), {3, 3});
    assert(AB.num_elements() == ab.size());

    for (int i = 0; i != C.size(); ++i, cout << '\n')
      for (int j = 0; j != std::get<1>(C.sizes()); ++j)
        cout << C[i][j] << ' ';
    cout << '\n';

    assert(C == AB);


    using ma::N;
    ma::product(N(A), N(B), C); // same as ma::product(A, B, C);
    assert(C == AB);

    using ma::T;

    ma::product(T(A), B, C);
    std::vector<double> atb = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
    boost::multi::array_ref<double, 2> AtB(atb.data(), {3, 3});
    assert(C == AtB);

    ma::product(A, T(B), C);
    std::vector<double> abt = {14., 14., 10., 92., 92., 110., 112., 121., 141.};
    boost::multi::array_ref<double, 2> ABt(abt.data(), {3, 3});
    assert(C == ABt);

    ma::product(T(A), T(B), C);
    std::vector<double> atbt = {44., 44., 58., 74., 65., 107., 94., 94., 138.};
    boost::multi::array_ref<double, 2> AtBt(atbt.data(), {3, 3});
    assert(C == AtBt);
  }

  cout << "test ended" << std::endl;
}
#endif

#endif
