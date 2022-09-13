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

#ifndef MA_LAPACK_HPP
#define MA_LAPACK_HPP

#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Numerics/detail/lapack.hpp"
#include "multi/array_ref.hpp"
#include <cassert>

namespace ma
{
template<class MultiArray2D>
int getrf_optimal_workspace_size(MultiArray2D&& A)
{
  assert(A.stride(0) > 0);
  assert(A.stride(1) == 1);

  int res;
  getrf_bufferSize(std::get<1>(A.sizes()), std::get<0>(A.sizes()), pointer_dispatch(A.origin()), A.stride(0), res);
  return res;
}

template<class MultiArray2D, class Array1D, class Buffer>
MultiArray2D&& getrf(MultiArray2D&& m, Array1D& pivot, Buffer&& WORK)
{
  assert(m.stride(0) >= std::max(std::size_t(1), std::size_t(std::get<1>(m.sizes()))));
  assert(m.stride(1) == 1);
  assert(pivot.size() >= std::min(std::get<1>(m.sizes()), std::get<0>(m.sizes()) + 1));

  int status = -1;
  getrf(std::get<1>(m.sizes()), std::get<0>(m.sizes()), pointer_dispatch(m.origin()), m.stride(0), pointer_dispatch(pivot.data()), status,
        pointer_dispatch(WORK.data()));
  //assert(status==0);
  return std::forward<MultiArray2D>(m);
}

template<class MultiArray2D>
int getri_optimal_workspace_size(MultiArray2D&& A)
{
  assert(A.stride(1) == 1);
  assert(std::get<0>(A.sizes()) ==std::get<1>(A.sizes()));
  int lwork = -1;
  getri_bufferSize(A.size(), pointer_dispatch(A.origin()), A.stride(), lwork);
  return lwork;
}

template<class MultiArray2D, class MultiArray1D, class Buffer>
MultiArray2D&& getri(MultiArray2D&& A, MultiArray1D const& IPIV, Buffer&& WORK)
{
  //	assert(A.stride(0) > std::max(std::size_t(1), A.size(1)));
  assert(A.stride(1) == 1);
  assert(IPIV.size() >= size_t(A.size()));
  assert(WORK.size() >= std::max(std::size_t(1), size_t(A.size())));

  int status = -1;
  getri(A.size(), pointer_dispatch(A.origin()), A.stride(), pointer_dispatch(IPIV.data()),
        pointer_dispatch(WORK.data()), WORK.size(), status);
  assert(status == 0);
  return std::forward<MultiArray2D>(A);
}

template<class MultiArray2D>
int geqrf_optimal_workspace_size(MultiArray2D&& A)
{
  assert(A.stride(0) > 0);
  assert(A.stride(1) == 1);

  int res;
  geqrf_bufferSize(std::get<1>(A.sizes()), std::get<0>(A.sizes()), pointer_dispatch(A.origin()), A.stride(0), res);
  return res;
}

template<class MultiArray2D, class Array1D, class Buffer>
MultiArray2D&& geqrf(MultiArray2D&& A, Array1D&& TAU, Buffer&& WORK)
{
  // why was this here???
  //assert(A.stride(0) > std::max(std::size_t(1), A.size(0)));
  assert(A.stride(1) == 1);
  assert(TAU.stride(0) == 1);
  assert(TAU.size() >= std::max(std::size_t(1), size_t(std::min(std::get<0>(A.sizes()), std::get<1>(A.sizes())))));
  assert(WORK.size() >= std::max(std::size_t(1), size_t(A.size())));

  int status = -1;
  geqrf(std::get<1>(A.sizes()), std::get<0>(A.sizes()), pointer_dispatch(A.origin()), A.stride(0), pointer_dispatch(TAU.origin()),
        pointer_dispatch(WORK.data()), WORK.size(), status);
  assert(status == 0);
  return std::forward<MultiArray2D>(A);
}

template<class MultiArray2D>
int gelqf_optimal_workspace_size(MultiArray2D&& A)
{
  assert(A.stride(0) > 0);
  assert(A.stride(1) == 1);

  int res;
  gelqf_bufferSize(std::get<1>(A.sizes()), std::get<0>(A.sizes()), pointer_dispatch(A.origin()), A.stride(0), res);
  return res;
}

template<class MultiArray2D, class Array1D, class Buffer>
MultiArray2D&& gelqf(MultiArray2D&& A, Array1D&& TAU, Buffer&& WORK)
{
  assert(A.stride(1) > 0);
  assert(A.stride(1) == 1);
  assert(TAU.stride(0) == 1);
  assert(TAU.size() >= std::max(std::size_t(1), size_t(std::min(std::get<0>(A.sizes()), std::get<1>(A.sizes())))));
  assert(WORK.size() >= std::max(std::size_t(1), size_t(std::get<1>(A.sizes()))));

  int status = -1;
  gelqf(std::get<1>(A.sizes()), std::get<0>(A.sizes()), pointer_dispatch(A.origin()), A.stride(0), pointer_dispatch(TAU.data()),
        pointer_dispatch(WORK.data()), WORK.size(), status);
  assert(status == 0);
  return std::forward<MultiArray2D>(A);
}


template<class MultiArray2D>
int gqr_optimal_workspace_size(MultiArray2D&& A)
{
  assert(A.stride(0) > 0);
  assert(A.stride(1) == 1);

  int res;
  gqr_bufferSize(std::get<1>(A.sizes()), std::get<0>(A.sizes()), std::max(std::size_t(1), size_t(std::min(std::get<0>(A.sizes()), std::get<1>(A.sizes())))),
                 pointer_dispatch(A.origin()), A.stride(0), res);
  return res;
}

template<class MultiArray2D, class Array1D, class Buffer>
MultiArray2D&& gqr(MultiArray2D&& A, Array1D&& TAU, Buffer&& WORK)
{
  assert(A.stride(1) == 1);
  assert(TAU.stride(0) == 1);
  assert(TAU.size() >= std::max(std::size_t(1), size_t(std::min(std::get<0>(A.sizes()), std::get<1>(A.sizes())))));
  assert(WORK.size() >= std::max(std::size_t(1), size_t(A.size())));

  int status = -1;
  gqr(std::get<1>(A.sizes()), std::get<0>(A.sizes()), std::max(std::size_t(1), size_t(std::min(std::get<0>(A.sizes()), std::get<1>(A.sizes())))),
      pointer_dispatch(A.origin()), A.stride(0), pointer_dispatch(TAU.origin()), pointer_dispatch(WORK.data()),
      WORK.size(), status);
  assert(status == 0);
  return std::forward<MultiArray2D>(A);
}

template<class MultiArray2D>
int glq_optimal_workspace_size(MultiArray2D&& A)
{
  assert(A.stride(0) > 0);
  assert(A.stride(1) == 1);

  int res;
  glq_bufferSize(std::get<1>(A.sizes()), std::get<0>(A.sizes()), std::max(std::size_t(1), size_t(std::min(std::get<0>(A.sizes()), std::get<1>(A.sizes())))),
                 pointer_dispatch(A.origin()), A.stride(0), res);
  return res;
}

template<class MultiArray2D, class Array1D, class Buffer>
MultiArray2D&& glq(MultiArray2D&& A, Array1D&& TAU, Buffer&& WORK)
{
  assert(A.stride(1) == 1);
  assert(TAU.stride(0) == 1);
  assert(TAU.size() >= std::max(std::size_t(1), size_t(std::min(std::get<0>(A.sizes()), std::get<1>(A.sizes())))));
  assert(WORK.size() >= std::max(std::size_t(1), size_t(std::get<1>(A.sizes()))));

  int status = -1;
  glq(std::get<1>(A.sizes()), std::get<0>(A.sizes()), std::max(std::size_t(1), size_t(std::min(std::get<0>(A.sizes()), std::get<1>(A.sizes())))),
      pointer_dispatch(A.origin()), A.stride(0), pointer_dispatch(TAU.data()), pointer_dispatch(WORK.data()),
      WORK.size(), status);
  assert(status == 0);
  return std::forward<MultiArray2D>(A);
}

template<class MultiArray2D, typename = typename std::enable_if_t<MultiArray2D::dimensionality == 2>>
MultiArray2D&& potrf(MultiArray2D&& A)
{
  assert(std::get<0>(A.sizes()) == std::get<1>(A.sizes()));
  int INFO;
  potrf('U', A.size(), pointer_dispatch(A.origin()), A.stride(0), INFO);
  if (INFO != 0)
    throw std::runtime_error(" error in ma::potrf: Error code != 0");
}

template<class MultiArray2D>
int gesvd_optimal_workspace_size(MultiArray2D&& A)
{
  assert(A.stride(0) > 0);
  assert(A.stride(1) == 1);

  int res;
  gesvd_bufferSize(std::get<1>(A.sizes()), std::get<0>(A.sizes()), pointer_dispatch(A.origin()), res);
  return res;
}

template<class MultiArray2D, class Array1D, class MultiArray2DU, class MultiArray2DV, class Buffer, class RBuffer>
MultiArray2D&& gesvd(char jobU,
                     char jobVT,
                     MultiArray2D&& A,
                     Array1D&& S,
                     MultiArray2DU&& U,
                     MultiArray2DV&& VT,
                     Buffer&& WORK,
                     RBuffer&& RWORK)
{
  assert(A.stride(1) > 0);
  assert(A.stride(1) == 1);

  // in C: A = U * S * VT
  // in F: At = (U * S * VT)t = VTt * S * Ut
  // so I need to switch U <--> VT when calling fortran interface
  int status = -1;
  gesvd(jobVT, jobU, std::get<1>(A.sizes()), std::get<0>(A.sizes()), pointer_dispatch(A.origin()), A.stride(0), pointer_dispatch(S.origin()),
        pointer_dispatch(VT.origin()), VT.stride(0), // !!!
        pointer_dispatch(U.origin()), U.stride(0),   // !!!
        pointer_dispatch(WORK.data()), WORK.size(), pointer_dispatch(RWORK.origin()), status);
  assert(status == 0);
  return std::forward<MultiArray2D>(A);
}

template<class MultiArray1D,
         class MultiArray2D,
         typename = typename std::enable_if_t<MultiArray1D::dimensionality == 1>,
         typename = typename std::enable_if_t<MultiArray2D::dimensionality == 2>>
std::pair<MultiArray1D, MultiArray2D> symEig(MultiArray2D const& A)
{
  using eigSys     = std::pair<MultiArray1D, MultiArray2D>;
  using Type       = typename MultiArray2D::element;
  using RealType   = typename qmcplusplus::afqmc::remove_complex<Type>::value_type;
  using extensions = typename boost::multi::layout_t<1u>::extensions_type;
  assert(A.size() == std::get<1>(A.sizes()));
  assert(A.stride(1) == 1);
  assert(A.size() > 0);
  int N   = A.size();
  int LDA = A.stride();

  MultiArray1D eigVal(extensions{N});
  MultiArray2D eigVec({N, N});
  MultiArray2D A_({N, N});
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      A_[i][j] = ma::conj(A[i][j]);
  char JOBZ('V');
  char RANGE('A');
  char UPLO('U');
  RealType VL     = 0;
  RealType VU     = 0;
  int IL          = 0;
  int IU          = 0;
  RealType ABSTOL = 0; //DLAMCH( 'Safe minimum' );
  int M;               // output: total number of eigenvalues found
  std::vector<int> ISUPPZ(2 * N);
  std::vector<Type> WORK(1); // set with workspace query
  int LWORK = -1;
  std::vector<RealType> RWORK(1); // set with workspace query
  int LRWORK = -1;
  std::vector<int> IWORK(1);
  int LIWORK = -1;
  int INFO;

  hevr(JOBZ, RANGE, UPLO, N, pointer_dispatch(A_.origin()), LDA, VL, VU, IL, IU, ABSTOL, M,
       pointer_dispatch(eigVal.origin()), pointer_dispatch(eigVec.origin()), N, pointer_dispatch(ISUPPZ.data()),
       pointer_dispatch(WORK.data()), LWORK, pointer_dispatch(RWORK.data()), LRWORK, pointer_dispatch(IWORK.data()),
       LIWORK, INFO);

  LWORK = int(real(WORK[0]));
  WORK.resize(LWORK);
  LRWORK = int(RWORK[0]);
  RWORK.resize(LRWORK);
  LIWORK = int(IWORK[0]);
  IWORK.resize(LIWORK);

  hevr(JOBZ, RANGE, UPLO, N, pointer_dispatch(A_.origin()), LDA, VL, VU, IL, IU, ABSTOL, M,
       pointer_dispatch(eigVal.origin()), pointer_dispatch(eigVec.origin()), N, pointer_dispatch(ISUPPZ.data()),
       pointer_dispatch(WORK.data()), LWORK, pointer_dispatch(RWORK.data()), LRWORK, pointer_dispatch(IWORK.data()),
       LIWORK, INFO);
  if (INFO != 0)
    throw std::runtime_error(" error in ma::eig: Error code != 0");
  if (M != N)
    throw std::runtime_error(" error in ma::eig: Not enough eigenvalues");
  for (int i = 0; i < N; i++)
    for (int j = i + 1; j < N; j++)
      std::swap(eigVec[i][j], eigVec[j][i]);

  return std::pair<MultiArray1D, MultiArray2D>{eigVal, eigVec};
}

// Careful!!!!
// This routine modifies the original matrix.
template<class MultiArray1D,
         class MultiArray2D,
         class MultiArray2DA,
         typename = typename std::enable_if_t<MultiArray1D::dimensionality == 1>,
         typename = typename std::enable_if_t<MultiArray2DA::dimensionality == 2>,
         typename = typename std::enable_if_t<MultiArray2D::dimensionality == 2>>
std::pair<MultiArray1D, MultiArray2D> symEigSelect(MultiArray2DA& A, int neig)
{
  using eigSys = std::pair<MultiArray1D, MultiArray2D>;
  using Type   = typename MultiArray2D::element;
  using TypeA  = typename MultiArray2DA::element;
  static_assert(std::is_same<Type, TypeA>::value, "Wrong types.");
  using RealType   = typename qmcplusplus::afqmc::remove_complex<Type>::value_type;
  using extensions = typename boost::multi::layout_t<1u>::extensions_type;
  assert(std::get<0>(A.sizes()) == std::get<1>(A.sizes()));
  assert(A.stride(1) == 1);
  assert(std::get<0>(A.sizes()) > 0);
  int N   = std::get<0>(A.sizes());
  int LDA = A.stride(0);

  MultiArray1D eigVal(extensions{neig});
  MultiArray2D eigVec({neig, N});
  // Transpose A to avoid using more memory
  for (int i = 0; i < N; i++)
    for (int j = i + 1; j < N; j++)
    {
      using std::swap;
      swap(A[i][j], A[j][i]);
    }

  char JOBZ('V');
  char RANGE('I');
  char UPLO('U');
  RealType VL     = 0;
  RealType VU     = 0;
  int IL          = 1;
  int IU          = neig;
  RealType ABSTOL = 0; //DLAMCH( 'Safe minimum' );
  int M;               // output: total number of eigenvalues found
  std::vector<int> ISUPPZ(2 * N);
  std::vector<Type> WORK(1); // set with workspace query
  int LWORK = -1;
  std::vector<RealType> RWORK(1); // set with workspace query
  int LRWORK = -1;
  std::vector<int> IWORK(1);
  int LIWORK = -1;
  int INFO;

  hevr(JOBZ, RANGE, UPLO, N, pointer_dispatch(A.origin()), LDA, VL, VU, IL, IU, ABSTOL, M,
       pointer_dispatch(eigVal.origin()), pointer_dispatch(eigVec.origin()), N, pointer_dispatch(ISUPPZ.data()),
       pointer_dispatch(WORK.data()), LWORK, pointer_dispatch(RWORK.data()), LRWORK, pointer_dispatch(IWORK.data()),
       LIWORK, INFO);

  LWORK = int(real(WORK[0]));
  WORK.resize(LWORK);
  LRWORK = int(RWORK[0]);
  RWORK.resize(LRWORK);
  LIWORK = int(IWORK[0]);
  IWORK.resize(LIWORK);

  hevr(JOBZ, RANGE, UPLO, N, pointer_dispatch(A.origin()), LDA, VL, VU, IL, IU, ABSTOL, M,
       pointer_dispatch(eigVal.origin()), pointer_dispatch(eigVec.origin()), N, pointer_dispatch(ISUPPZ.data()),
       pointer_dispatch(WORK.data()), LWORK, pointer_dispatch(RWORK.data()), LRWORK, pointer_dispatch(IWORK.data()),
       LIWORK, INFO);
  if (INFO != 0)
    throw std::runtime_error(" error in ma::eig: Error code != 0");
  if (M != neig)
    throw std::runtime_error(" error in ma::eig: Not enough eigenvalues");

  return std::pair<MultiArray1D, MultiArray2D>{eigVal, eigVec};
}

// Careful!!!!
// This routine modifies the original matrix.
template<class MultiArray1D,
         class MultiArray2D,
         class MultiArray2DA,
         class MultiArray2DB,
         typename = typename std::enable_if_t<MultiArray1D::dimensionality == 1>,
         typename = typename std::enable_if_t<MultiArray2DA::dimensionality == 2>,
         typename = typename std::enable_if_t<MultiArray2DB::dimensionality == 2>,
         typename = typename std::enable_if_t<MultiArray2D::dimensionality == 2>>
std::pair<MultiArray1D, MultiArray2D> genEigSelect(MultiArray2DA& A, MultiArray2DB& S, int neig, int itype = 1)
{
  using eigSys = std::pair<MultiArray1D, MultiArray2D>;
  using Type   = typename MultiArray2D::element;
  using TypeA  = typename MultiArray2DA::element;
  using TypeB  = typename MultiArray2DB::element;
  static_assert(std::is_same<Type, TypeA>::value, "Wrong types.");
  static_assert(std::is_same<TypeA, TypeB>::value, "Wrong types.");
  using RealType   = typename qmcplusplus::afqmc::remove_complex<Type>::value_type;
  using extensions = typename boost::multi::layout_t<1u>::extensions_type;
  assert(std::get<0>(A.sizes()) == std::get<1>(A.sizes()));
  assert(std::get<0>(A.sizes()) == std::get<0>(S.sizes()));
  assert(std::get<0>(S.sizes()) == std::get<1>(S.sizes()));
  assert(A.stride(1) == 1);
  assert(std::get<0>(A.sizes()) > 0);
  assert(S.stride(1) == 1);
  assert(std::get<0>(S.sizes()) > 0);
  int N   = std::get<0>(A.sizes());
  int LDA = A.stride(0);
  int LDS = S.stride(0);

  MultiArray1D eigVal(extensions{neig});
  MultiArray2D eigVec({neig, N});
  // Transpose A to avoid using more memory
  for (int i = 0; i < N; i++)
    for (int j = i + 1; j < N; j++)
    {
      using std::swap;
      swap(A[i][j], A[j][i]);
      swap(S[i][j], S[j][i]);
    }

  char JOBZ('V');
  char RANGE('I');
  char UPLO('U');
  RealType VL     = 0;
  RealType VU     = 0;
  int IL          = 1;
  int IU          = neig;
  RealType ABSTOL = 0;       //DLAMCH( 'Safe minimum' );
  int M;                     // output: total number of eigenvalues found
  std::vector<Type> WORK(1); // set with workspace query
  int LWORK = -1;
  std::vector<RealType> RWORK(7 * N); // set with workspace query
  std::vector<int> IWORK(5 * N);
  std::vector<int> IFAIL(N);
  int INFO;

  gvx(itype, JOBZ, RANGE, UPLO, N, pointer_dispatch(A.origin()), LDA, pointer_dispatch(S.origin()), LDS, VL, VU, IL, IU,
      ABSTOL, M, pointer_dispatch(eigVal.origin()), pointer_dispatch(eigVec.origin()), N, pointer_dispatch(WORK.data()),
      LWORK, pointer_dispatch(RWORK.data()), pointer_dispatch(IWORK.data()), IFAIL.data(), INFO);

  LWORK = int(real(WORK[0]));
  WORK.resize(LWORK);

  gvx(itype, JOBZ, RANGE, UPLO, N, pointer_dispatch(A.origin()), LDA, pointer_dispatch(S.origin()), LDS, VL, VU, IL, IU,
      ABSTOL, M, pointer_dispatch(eigVal.origin()), pointer_dispatch(eigVec.origin()), N, pointer_dispatch(WORK.data()),
      LWORK, pointer_dispatch(RWORK.data()), pointer_dispatch(IWORK.data()), IFAIL.data(), INFO);
  if (INFO != 0)
    throw std::runtime_error(" error in ma::eig: Error code != 0");
  if (M != neig)
    throw std::runtime_error(" error in ma::eig: Not enough eigenvalues");

  return std::pair<MultiArray1D, MultiArray2D>{eigVal, eigVec};
}

} // namespace ma

#endif
