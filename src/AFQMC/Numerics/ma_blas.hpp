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

#ifndef MA_BLAS_HPP
#define MA_BLAS_HPP

#include "AFQMC/Numerics/detail/blas.hpp"
#include <utility> //std::enable_if
#include <cassert>
#include <iostream>

#include "AFQMC/Numerics/ma_blas_extensions.hpp"

namespace ma
{
template<class MultiArray1DX,
         class MultiArray1DY,
         typename = typename std::enable_if_t<std::decay<MultiArray1DX>::type::dimensionality == 1>,
         typename = typename std::enable_if_t<std::decay<MultiArray1DY>::type::dimensionality == 1>>
MultiArray1DY&& copy(MultiArray1DX&& x, MultiArray1DY&& y)
{
  assert(x.num_elements() == y.num_elements());
  copy(x.size(), pointer_dispatch(x.origin()), x.stride(), pointer_dispatch(y.origin()), y.stride());
  return std::forward<MultiArray1DY>(y);
}

template<class MultiArray2DX,
         class MultiArray2DY,
         typename = typename std::enable_if_t<std::decay<MultiArray2DX>::type::dimensionality == 2>,
         typename = typename std::enable_if_t<std::decay<MultiArray2DY>::type::dimensionality == 2>,
         typename = void>
MultiArray2DY&& copy(MultiArray2DX&& x, MultiArray2DY&& y)
{
  assert(x.stride(1) == 1);
  assert(y.stride(1) == 1);
  assert(x.size() == y.size());
  assert(std::get<1>(x.sizes()) == std::get<1>(y.sizes()));
  if ((x.stride() == std::get<1>(x.sizes())) && (y.stride() == std::get<1>(y.sizes())))
  {
    copy(x.num_elements(), pointer_dispatch(x.origin()), 1, pointer_dispatch(y.origin()), 1);
  }
  else
  {
    copy2D(x.size(), std::get<1>(x.sizes()), pointer_dispatch(x.origin()), x.stride(), pointer_dispatch(y.origin()), y.stride());
  }
  return std::forward<MultiArray2DY>(y);
}

template<class MultiArrayNDX,
         class MultiArrayNDY,
         typename = typename std::enable_if_t<(std::decay<MultiArrayNDX>::type::dimensionality > 2)>,
         typename = typename std::enable_if_t<(std::decay<MultiArrayNDY>::type::dimensionality > 2)>,
         typename = void,
         typename = void>
MultiArrayNDY&& copy(MultiArrayNDX&& x, MultiArrayNDY&& y)
{
#ifndef NDEBUG
  // only on contiguous arrays
//  long sz(x.size());
//  for (int i = 1; i < int(std::decay<MultiArrayNDX>::type::dimensionality); ++i)
//    sz *= x.size(i);
//  assert(x.num_elements() == sz);
  assert(std::get<std::decay<MultiArrayNDX>::type::dimensionality - 1>(x.strides()) == 1);
//  sz = y.size();
//  for (int i = 1; i < int(std::decay<MultiArrayNDY>::type::dimensionality); ++i)
//   sz *= y.size(i);
//  assert(y.num_elements() == sz);
  assert(std::get<std::decay<MultiArrayNDY>::type::dimensionality - 1>(y.strides()) == 1);
  assert(x.num_elements() == y.num_elements());
#endif
  copy(x.num_elements(), pointer_dispatch(x.origin()), 1, pointer_dispatch(y.origin()), 1);
  return std::forward<MultiArrayNDY>(y);
}

template<class MultiArray1Dx,
         class MultiArray1Dy,
         typename = typename std::enable_if<std::decay<MultiArray1Dx>::type::dimensionality == 1>::type,
         typename = typename std::enable_if<std::decay<MultiArray1Dy>::type::dimensionality == 1>::type>
//auto
typename std::decay<MultiArray1Dx>::type::element dot(MultiArray1Dx&& x, MultiArray1Dy&& y)
{
  assert(x.size() == y.size());
  return dot(x.size(), pointer_dispatch(x.origin()), x.stride(), pointer_dispatch(y.origin()), y.stride());
}

template<class MultiArray2Dx,
         class MultiArray2Dy,
         typename = typename std::enable_if<std::decay<MultiArray2Dx>::type::dimensionality == 2>::type,
         typename = typename std::enable_if<std::decay<MultiArray2Dy>::type::dimensionality == 2>::type,
         typename = void>
typename std::decay<MultiArray2Dx>::type::element dot(MultiArray2Dx&& x, MultiArray2Dy&& y)
{
  assert(x.stride() == std::get<1>(x.sizes())); // only on contiguous arrays
  assert(x.stride(1) == 1);         // only on contiguous arrays
  assert(y.stride() == std::get<1>(y.sizes())); // only on contiguous arrays
  assert(y.stride(1) == 1);         // only on contiguous arrays
  assert(x.num_elements() == y.num_elements());
  return dot(x.num_elements(), pointer_dispatch(x.origin()), 1, pointer_dispatch(y.origin()), 1);
}

template<class T,
         class MultiArray1D,
         typename = typename std::enable_if<std::decay<MultiArray1D>::type::dimensionality == 1>::type>
MultiArray1D&& scal(T a, MultiArray1D&& x)
{
  scal(x.size(), a, pointer_dispatch(x.origin()), x.stride());
  return std::forward<MultiArray1D>(x);
}

template<class T,
         class MultiArrayND,
         typename = typename std::enable_if<(std::decay<MultiArrayND>::type::dimensionality > 1)>::type,
         typename = void // TODO change to use dispatch
         >
MultiArrayND&& scal(T a, MultiArrayND&& x)
{
#ifndef NDEBUG
//  long sz(x.size());
//  for (int i = 1; i < int(std::decay<MultiArrayND>::type::dimensionality); ++i)
//    sz *= x.size(i);
//  assert(x.num_elements() == sz);
  assert(std::get<std::decay<MultiArrayND>::type::dimensionality - 1>(x.strides()) == 1); // only on contiguous arrays
#endif
  scal(x.num_elements(), a, pointer_dispatch(x.origin()), 1);
  return std::forward<MultiArrayND>(x);
}
/*
template<class T,
        class MultiArray3D,
        typename = typename std::enable_if<std::decay<MultiArray3D>::type::dimensionality == 3>::type,
        typename = void, // TODO change to use dispatch 
        typename = void // TODO change to use dispatch 
    >
MultiArray3D scal(T a, MultiArray3D&& x){
        assert( x.stride(0) == x.size(1)*x.size(2) ); // only on contiguous arrays 
        assert( x.stride(1) == x.size(2) ); // only on contiguous arrays 
        assert( x.stride(2) == 1 );            // only on contiguous arrays 
        scal(x.num_elements(), a, pointer_dispatch(x.origin()), 1);
        return std::forward<MultiArray3D>(x);
}
*/
template<class T, class MultiArray1D>
auto operator*=(MultiArray1D&& x, T a) -> decltype(scal(a, std::forward<MultiArray1D>(x)))
{
  return scal(a, std::forward<MultiArray1D>(x));
}

template<class T,
         class MultiArray1DA,
         class MultiArray1DB,
         typename = typename std::enable_if<std::decay<MultiArray1DA>::type::dimensionality == 1 and
                                            std::decay<MultiArray1DB>::type::dimensionality == 1>::type>
MultiArray1DB&& axpy(T x, MultiArray1DA const& a, MultiArray1DB&& b)
{
  assert(a.size() == b.size());
  axpy(a.size(), x, pointer_dispatch(a.origin()), a.stride(), pointer_dispatch(b.origin()), b.stride());
  return std::forward<MultiArray1DB>(b);
}

template<class T,
         class MultiArray2DA,
         class MultiArray2DB,
         typename = typename std::enable_if<std::decay<MultiArray2DA>::type::dimensionality == 2 and
                                            std::decay<MultiArray2DB>::type::dimensionality == 2>::type,
         typename = void // TODO change to use dispatch
         >
MultiArray2DB&& axpy(T x, MultiArray2DA const& a, MultiArray2DB&& b)
{
  assert(a.num_elements() == b.num_elements());
  assert(a.stride() == std::get<1>(a.sizes())); // only on contiguous arrays
  assert(a.stride(1) == 1);         // only on contiguous arrays
  assert(b.stride() == std::get<1>(b.sizes())); // only on contiguous arrays
  assert(b.stride(1) == 1);         // only on contiguous arrays
  axpy(a.num_elements(), x, pointer_dispatch(a.origin()), 1, pointer_dispatch(b.origin()), 1);
  return std::forward<MultiArray2DB>(b);
}

template<
    char IN,
    class T,
    class MultiArray2DA,
    class MultiArray1DX,
    class MultiArray1DY,
    typename = typename std::enable_if<MultiArray2DA::dimensionality == 2 and MultiArray1DX::dimensionality == 1 and
                                       std::decay<MultiArray1DY>::type::dimensionality == 1>::type>
MultiArray1DY&& gemv(T alpha, MultiArray2DA const& A, MultiArray1DX const& x, T beta, MultiArray1DY&& y)
{
  assert((IN == 'N') || (IN == 'T') || (IN == 'C'));
  if (IN == 'T' or IN == 'C')
    assert(x.size() == std::get<1>(A.sizes()) and y.size() == A.size());
  else if (IN == 'N')
    assert(x.size() == A.size() and y.size() == std::get<1>(A.sizes()));
  assert(A.stride(1) == 1); // gemv is not implemented for arrays with non-leading stride != 1
  int M = std::get<1>(A.sizes());
  int N = A.size();
  gemv(IN, M, N, alpha, pointer_dispatch(A.origin()), A.stride(), pointer_dispatch(x.origin()), x.stride(), beta,
       pointer_dispatch(y.origin()), y.stride());
  return std::forward<MultiArray1DY>(y);
} //y := alpha*A*x + beta*y,

template<char IN, class MultiArray2DA, class MultiArray1DX, class MultiArray1DY>
MultiArray1DY&& gemv(MultiArray2DA const& A, MultiArray1DX const& x, MultiArray1DY&& y)
{
  return gemv<IN>(1., A, x, 0., std::forward<MultiArray1DY>(y));
} //y := alpha*A*x

//	gemm<'T', 'T'>(1., A, B, 0., C); // C = T(A*B) = T(B)*T(A) or T(C) = A*B
//	gemm<'N', 'N'>(1., A, B, 0., C); // C = B*A = T(T(A)*T(B)) or T(C) = T(A)*T(B)
//	gemm<'T', 'N'>(1., A, B, 0., C); // C = T(A*T(B)) = B*T(A) or T(C) = A*T(B)
//	gemm<'N', 'T'>(1., A, B, 0., C); // C =  T(T(A)*B) = T(B)*A or T(C) = T(A)*B

template<
    char TA,
    char TB,
    class T,
    class MultiArray2DA,
    class MultiArray2DB,
    class MultiArray2DC,
    typename = typename std::enable_if<MultiArray2DA::dimensionality == 2 and MultiArray2DB::dimensionality == 2 and
                                       std::decay<MultiArray2DC>::type::dimensionality == 2>::type>
MultiArray2DC&& gemm(T alpha, MultiArray2DA const& a, MultiArray2DB const& b, T beta, MultiArray2DC&& c)
{
  assert(a.stride(1) == 1);
  assert(b.stride(1) == 1);
  assert(c.stride(1) == 1);
  assert((TA == 'N') || (TA == 'T') || (TA == 'C'));
  assert((TB == 'N') || (TB == 'T') || (TB == 'C'));
  int M = -1;
  int N = -1;
  int K = -1;
  if (TA == 'N' and TB == 'N')
  {
    M = std::get<1>(a.sizes());
    N = b.size();
    K = a.size();
    assert(a.size() == std::get<1>(b.sizes()) and c.size() == b.size() and std::get<1>(c.sizes()) == std::get<1>(a.sizes()));
  }
  if ((TA == 'T' or TA == 'C') and (TB == 'T' or TB == 'C'))
  {
    M = a.size();
    N = std::get<1>(b.sizes());
    K = std::get<1>(a.sizes());
    assert(std::get<1>(a.sizes()) == b.size() and c.size() == std::get<1>(b.sizes()) and std::get<1>(c.sizes()) == a.size());
  }
  if ((TA == 'T' or TA == 'C') and TB == 'N')
  {
    M = a.size();
    N = b.size();
    K = std::get<1>(a.sizes());
    assert(std::get<1>(a.sizes()) == std::get<1>(b.sizes()) and c.size() == b.size() and std::get<1>(c.sizes()) == a.size());
  }
  if (TA == 'N' and (TB == 'T' or TB == 'C'))
  {
    M = std::get<1>(a.sizes());
    N = std::get<1>(b.sizes());
    K = a.size();
    assert(a.size() == b.size() and c.size() == std::get<1>(b.sizes()) and std::get<1>(c.sizes()) == std::get<1>(a.sizes()));
  }
  gemm(TA, TB, M, N, K, alpha, pointer_dispatch(a.origin()), a.stride(), pointer_dispatch(b.origin()), b.stride(),
       beta, pointer_dispatch(c.origin()), c.stride());
  return std::forward<MultiArray2DC>(c);
}

// Expect: A[nbatch][nrow][ncol]
template<
    char TA,
    char TB,
    class T,
    class MultiArray3DA,
    class MultiArray3DB,
    class MultiArray3DC,
    typename = typename std::enable_if<MultiArray3DA::dimensionality == 3 and MultiArray3DB::dimensionality == 3 and
                                       std::decay<MultiArray3DC>::type::dimensionality == 3>::type>
MultiArray3DC&& gemmStridedBatched(T alpha, MultiArray3DA const& a, MultiArray3DB const& b, T beta, MultiArray3DC&& c)
{
  assert(a.stride(2) == 1);
  assert(b.stride(2) == 1);
  assert(c.stride(2) == 1);
  assert(a.size() == b.size());
  assert(a.size() == c.size());
  assert((TA == 'N') || (TA == 'T') || (TA == 'C'));
  assert((TB == 'N') || (TB == 'T') || (TB == 'C'));
  int M = -1;
  int N = -1;
  int K = -1;
  if (TA == 'N' and TB == 'N')
  {
    M = std::get<2>(a.sizes());
    N = std::get<1>(b.sizes());
    K = std::get<1>(a.sizes());
    assert(std::get<1>(a.sizes()) == std::get<2>(b.sizes()) and std::get<1>(c.sizes()) == std::get<1>(b.sizes()) and std::get<2>(c.sizes()) == std::get<2>(a.sizes()));
  }
  if ((TA == 'T' or TA == 'C') and (TB == 'T' or TB == 'C'))
  {
    M = std::get<1>(a.sizes());
    N = std::get<2>(b.sizes());
    K = std::get<2>(a.sizes());
    assert(std::get<2>(a.sizes()) == std::get<1>(b.sizes()) and std::get<1>(c.sizes()) == std::get<2>(b.sizes()) and std::get<2>(c.sizes()) == std::get<1>(a.sizes()));
  }
  if ((TA == 'T' or TA == 'C') and TB == 'N')
  {
    M = std::get<1>(a.sizes());
    N = std::get<1>(b.sizes());
    K = std::get<2>(a.sizes());
    assert(std::get<2>(a.sizes()) == std::get<2>(b.sizes()) and std::get<1>(c.sizes()) == std::get<1>(b.sizes()) and std::get<2>(c.sizes()) == std::get<1>(a.sizes()));
  }
  if (TA == 'N' and (TB == 'T' or TB == 'C'))
  {
    M = std::get<2>(a.sizes());
    N = std::get<2>(b.sizes());
    K = std::get<1>(a.sizes());
    assert(std::get<1>(a.sizes()) == std::get<1>(b.sizes()) and std::get<1>(c.sizes()) == std::get<2>(b.sizes()) and std::get<2>(c.sizes()) == std::get<2>(a.sizes()));
  }
  gemmStridedBatched(TA, TB, M, N, K, alpha, pointer_dispatch(a.origin()), a.stride(1), a.stride(),
                     pointer_dispatch(b.origin()), b.stride(1), b.stride(), beta, pointer_dispatch(c.origin()),
                     c.stride(1), c.stride(), a.size());
  return std::forward<MultiArray3DC>(c);
}

template<char TA, char TB, class T, class MultiArray2DA, class MultiArray2DB, class MultiArray2DC>
MultiArray2DC&& gemm(MultiArray2DA const& a, MultiArray2DB const& b, MultiArray2DC&& c)
{
  return gemm(1., a, b, 0., std::forward<MultiArray2DC>(c));
}

template<
    char TA,
    char TB,
    class T,
    class MultiArray2DA,
    class MultiArray2DB,
    class MultiArray2DC,
    typename = typename std::enable_if<MultiArray2DA::dimensionality == 2 and MultiArray2DB::dimensionality == 2 and
                                       std::decay<MultiArray2DC>::type::dimensionality == 2>::type>
MultiArray2DC&& geam(T alpha, MultiArray2DA const& a, T beta, MultiArray2DB const& b, MultiArray2DC&& c)
{
  assert(a.stride(1) == 1);
  assert(b.stride(1) == 1);
  assert(c.stride(1) == 1);
  assert((TA == 'N') || (TA == 'T') || (TA == 'C'));
  assert((TB == 'N') || (TB == 'T') || (TB == 'C'));
  if (TA == 'N' and TB == 'N')
  {
    assert(a.size() == c.size() and std::get<1>(a.sizes()) == std::get<1>(c.sizes()));
    assert(b.size() == c.size() and std::get<1>(b.sizes()) == std::get<1>(c.sizes()));
  }
  if ((TA == 'T' or TA == 'C') and (TB == 'T' or TB == 'C'))
  {
    assert(std::get<1>(a.sizes()) == c.size() and a.size() == std::get<1>(c.sizes()));
    assert(std::get<1>(b.sizes()) == c.size() and b.size() == std::get<1>(c.sizes()));
  }
  if ((TA == 'T' or TA == 'C') and TB == 'N')
  {
    assert(std::get<1>(a.sizes()) == c.size() and a.size() == std::get<1>(c.sizes()));
    assert(b.size() == c.size() and std::get<1>(b.sizes()) == std::get<1>(c.sizes()));
  }
  if (TA == 'N' and (TB == 'T' or TB == 'C'))
  {
    assert(a.size() == c.size() and std::get<1>(a.sizes()) == std::get<1>(c.sizes()));
    assert(std::get<1>(b.sizes()) == c.size() and b.size() == std::get<1>(c.sizes()));
  }
  geam(TA, TB, std::get<1>(c.sizes()), c.size(), alpha, pointer_dispatch(a.origin()), a.stride(), beta,
       pointer_dispatch(b.origin()), b.stride(), pointer_dispatch(c.origin()), c.stride());
  return std::forward<MultiArray2DC>(c);
}

template<char TA,
         class T,
         class MultiArray2DA,
         class MultiArray2DC,
         typename = typename std::enable_if<MultiArray2DA::dimensionality == 2 and
                                            std::decay<MultiArray2DC>::type::dimensionality == 2>::type>
MultiArray2DC&& geam(T alpha, MultiArray2DA const& a, MultiArray2DC&& c)
{
  assert(a.stride(1) == 1);
  assert(c.stride(1) == 1);
  assert((TA == 'N') || (TA == 'T') || (TA == 'C'));
  if (TA == 'N')
  {
    assert(a.size() == c.size() and std::get<1>(a.sizes()) == std::get<1>(c.sizes()));
  }
  if ((TA == 'T' or TA == 'C'))
  {
    assert(std::get<1>(a.sizes()) == c.size() and a.size() == std::get<1>(c.sizes()));
  }
  geam(TA, TA, std::get<1>(c.sizes()), c.size(), alpha, pointer_dispatch(a.origin()), a.stride(), T(0),
       pointer_dispatch(a.origin()), a.stride(), pointer_dispatch(c.origin()), c.stride());
  return std::forward<MultiArray2DC>(c);
}

} // namespace ma

#endif
