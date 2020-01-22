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

#ifndef MA_BLAS_EXTENSIONS_HPP
#define MA_BLAS_EXTENSIONS_HPP

#include "AFQMC/Numerics/detail/blas.hpp"
#include<utility> //std::enable_if
#include<cassert>
#include<iostream>

namespace ma{

template<class T,
         class Q,
         class MultiArray1Dx,
         class MultiArray1Dy,
         class ptr,
         typename = typename std::enable_if<std::decay<MultiArray1Dx>::type::dimensionality == 1>::type,
         typename = typename std::enable_if<std::decay<MultiArray1Dy>::type::dimensionality == 1>::type
>
void adotpby(T const alpha, MultiArray1Dx const& x, MultiArray1Dy const& y, Q const beta, ptr res){
        assert(x.size() == y.size());
        adotpby(x.size(), alpha, pointer_dispatch(x.origin()), x.stride(0),
                                        pointer_dispatch(y.origin()), y.stride(0), beta, to_address(res));
}

template<class T,
         class MultiArray1Dx,
         class MultiArray1Dy,
         typename = typename std::enable_if<std::decay<MultiArray1Dx>::type::dimensionality == 1>::type,
         typename = typename std::enable_if<std::decay<MultiArray1Dy>::type::dimensionality == 1>::type
>
MultiArray1Dy
axty(T const alpha, MultiArray1Dx const& x, MultiArray1Dy && y){
        assert(x.size() == y.size());
        axty(x.size(), alpha, pointer_dispatch(x.origin()), x.stride(0), pointer_dispatch(y.origin()), y.stride(0));
        return y;
}

template<class T,
         class MultiArray2DA,
         class MultiArray2DB,
         typename = typename std::enable_if<std::decay<MultiArray2DA>::type::dimensionality == 2>::type,
         typename = typename std::enable_if<std::decay<MultiArray2DB>::type::dimensionality == 2>::type,
         typename = void
>
MultiArray2DB
axty(T const alpha, MultiArray2DA const& A, MultiArray2DB && B){
        assert(A.num_elements() == B.num_elements());
        assert(A.stride(1)==1);
        assert(A.stride(0)==A.size(1));
        assert(B.stride(1)==1);
        assert(B.stride(0)==B.size(1));
        axty(A.num_elements(), alpha, pointer_dispatch(A.origin()), 1, pointer_dispatch(B.origin()), 1);
        return B;
}

// on fortran ordering
// implements z[i][j] = beta * z[i][j] + alpha * conj(y[i][j]) * x[i]
template<class T,
         class MultiArray2DA,
         class MultiArray1D,
         class MultiArray2DB,
         typename = typename std::enable_if_t< (MultiArray2DA::dimensionality == 2) and
                                (MultiArray1D::dimensionality == 1) and
                                (std::decay<MultiArray2DB>::type::dimensionality == 2)>
>
MultiArray2DB
acAxpbB(T const alpha, MultiArray2DA const& A, MultiArray1D const& x, T const beta, MultiArray2DB && B){
        assert(A.num_elements() == B.num_elements());
        assert(A.size(0)==B.size(0));
        assert(A.size(1)==B.size(1));
        assert(A.size(1)==x.size(0));
        acAxpbB(A.size(1),A.size(0),alpha,pointer_dispatch(A.origin()),A.stride(0),
                pointer_dispatch(x.origin()),x.stride(0),beta,pointer_dispatch(B.origin()),B.stride(0));
        return B;
}

template<class T,
         class MultiArray2DA,
         class MultiArray1Dy,
         typename = typename std::enable_if<std::decay<MultiArray2DA>::type::dimensionality == 2>::type,
         typename = typename std::enable_if<std::decay<MultiArray1Dy>::type::dimensionality == 1>::type
>
MultiArray1Dy
adiagApy(T const alpha, MultiArray2DA const& A, MultiArray1Dy && y){
        assert(A.size(0) == A.size(1));
        assert(A.size(0) == y.size());
        adiagApy(y.size(), alpha, pointer_dispatch(A.origin()), A.stride(0), pointer_dispatch(y.origin()), y.stride(0));
        return y;
}

template<class MultiArray1D,
         typename = typename std::enable_if<std::decay<MultiArray1D>::type::dimensionality == 1>::type
>
auto
sum(MultiArray1D const& y){
        return sum(y.size(), pointer_dispatch(y.origin()), y.stride(0));
}

template<class MultiArray2D,
         typename = typename std::enable_if<std::decay<MultiArray2D>::type::dimensionality == 2>::type,
         typename = void
>
auto
sum(MultiArray2D const& A){
        assert(A.stride(1) == 1);
        // blas call assumes fortran ordering
        return sum(A.size(1), A.size(0), pointer_dispatch(A.origin()), A.stride(0));
}

template<class MultiArray3D,
         typename = typename std::enable_if<std::decay<MultiArray3D>::type::dimensionality == 3>::type,
         typename = void,
         typename = void
>
auto
sum(MultiArray3D const& A){
        // only arrays and array_refs for now
        assert(A.stride(0) == A.size(1)*A.size(2));
        assert(A.stride(1) == A.size(2));
        assert(A.stride(2) == 1);
        return sum(A.num_elements(), pointer_dispatch(A.origin()), 1);
}

template<class MultiArray4D,
         typename = typename std::enable_if<std::decay<MultiArray4D>::type::dimensionality == 4>::type,
         typename = void,
         typename = void,
         typename = void
>
auto
sum(MultiArray4D const& A){
        // only arrays and array_refs for now
        assert(A.stride(0) == A.size(1)*A.size(2)*A.size(3));
        assert(A.stride(1) == A.size(2)*A.size(3));
        assert(A.stride(2) == A.size(3));
        assert(A.stride(3) == 1);
        return sum(A.num_elements(), pointer_dispatch(A.origin()), 1);
}

template<class T, class MultiArray1D,
        typename = typename std::enable_if< std::decay<MultiArray1D>::type::dimensionality == 1 >
>
MultiArray1D setVector(T alpha, MultiArray1D&& a){
        set1D(a.size(0),  alpha, pointer_dispatch(a.origin()), a.stride(0) );
        return std::forward<MultiArray1D>(a);
}

template<class MultiArray1D,
        typename = std::enable_if_t< std::decay<MultiArray1D>::type::dimensionality == 1 >
>
void zero_complex_part(MultiArray1D&& a){
        zero_complex_part(a.num_elements(),pointer_dispatch(a.origin()));
}

template<class MultiArray2D,
        typename = std::enable_if_t< std::decay<MultiArray2D>::type::dimensionality == 2 >
        >
MultiArray2D set_identity(MultiArray2D&& m){
        set_identity(m.size(1),m.size(0),pointer_dispatch(m.origin()),m.stride(0));
        return std::forward<MultiArray2D>(m);
}

template<class MultiArray3D,
        typename = std::enable_if_t< std::decay<MultiArray3D>::type::dimensionality == 3 >,
        typename = void
        >
MultiArray3D set_identity(MultiArray3D&& m){
        set_identity_strided(m.size(0),m.stride(0),m.size(2),m.size(1),pointer_dispatch(m.origin()),m.stride(1));
        return std::forward<MultiArray3D>(m);
}

template<class T, class MultiArray2D,
        typename = typename std::enable_if< std::decay<MultiArray2D>::type::dimensionality == 2 >
>
MultiArray2D fill(MultiArray2D&& m, T const& value){
        using qmcplusplus::afqmc::fill2D;
        fill2D(m.size(0),m.size(1),pointer_dispatch(m.origin()),m.stride(0),value);
        return std::forward<MultiArray2D>(m);
}

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
  if(TA=='N' || TA=='Z') {
    if(M.size(0) != A.size(0) or M.size(1) != A.size(1)) 
      M.reextent({A.size(0),A.size(1)});
  } else if(TA=='T' || TA=='H') {
    if(M.size(0) != A.size(1) or M.size(1) != A.size(0))
      M.reextent({A.size(1),A.size(0)});
  } else {
    throw std::runtime_error(" Error: Unknown operation in Matrix2MA.\n");
  }
  using std::fill_n;
  fill_n(M.origin(),M.num_elements(),Type(0));
  auto pbegin = A.pointers_begin();
  auto pend = A.pointers_end();
  int_type p0(pbegin[0]);
  auto v0 = A.non_zero_values_data();
  auto c0 = A.non_zero_indices2_data();
//#ifdef ENABLE_CUDA
//qmcplusplus::app_log()<<" /**********************************\n";
//qmcplusplus::app_log()<<" Warning: write kernel in Matrix2MA. \n";
//qmcplusplus::app_log()<<" /**********************************\n";
//#endif
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
//#ifdef ENABLE_CUDA
//qmcplusplus::app_log()<<" /**********************************\n";
//qmcplusplus::app_log()<<" Warning: write kernel in Matrix2MAREF. \n";
//qmcplusplus::app_log()<<" /**********************************\n";
//#endif
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

/* Chooses rows of A based on occups vector and performs CSF2MA on subset of rows */
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

template<class MA,
         class MultiArray2D,
         typename = typename std::enable_if_t<(std::decay<MA>::type::dimensionality == 2)>,
         typename = typename std::enable_if_t<(MultiArray2D::dimensionality==2)>,
         typename = void
        >
void Matrix2MA(char TA, MA const& A, MultiArray2D& M)
{
  using Type1 = typename std::decay<MA>::type::element;
  using Type2 = typename MultiArray2D::element;
  assert(TA=='N' || TA=='H' || TA=='T' || TA=='Z');
  if(TA=='N' || TA=='Z') {
    if(M.size(0) != A.size(0) or M.size(1) != A.size(1)) 
      M.reextent({A.size(0),A.size(1)});
  } else if(TA=='T' || TA=='H') {
    if(M.size(0) != A.size(1) or M.size(1) != A.size(0)) 
      M.reextent({A.size(1),A.size(0)});
  } else {
    throw std::runtime_error(" Error: Unknown operation in Matrix2MA.\n");
  }

  using ptrA = std::remove_cv_t<typename MA::element_ptr>;
  using ptrM = std::remove_cv_t<typename MultiArray2D::element_ptr>;

  if(TA == 'H') TA = 'C';   
  if(TA=='Z') {
    for(int i=0; i<M.size(0); i++)
      for(int j=0; j<M.size(1); j++)
        M[i][j] = ma::conj(A[i][j]);
  } else if( not std::is_same<ptrA,ptrM>::value ) {
    if(TA=='N') {
      for(int i=0; i<M.size(0); i++)
        for(int j=0; j<M.size(1); j++)
          M[i][j] = A[i][j];
    } else if(TA=='T') {
      for(int i=0; i<M.size(0); i++)
        for(int j=0; j<M.size(1); j++)
          M[i][j] = A[j][i];
    } else if(TA=='C') {
      for(int i=0; i<M.size(0); i++)
        for(int j=0; j<M.size(1); j++)
          M[i][j] = ma::conj(A[j][i]);
    }
  } else {
    geam(   TA, TA, M.size(1), M.size(0),
        Type2(1.0), pointer_dispatch(A.origin()), A.stride(0),
        Type2(0.0), pointer_dispatch(A.origin()), A.stride(0),
        pointer_dispatch(M.origin()), M.stride(0)
    );
  }
}

template<class MA,
         class MultiArray2D,
         typename = typename std::enable_if_t<(std::decay<MA>::type::dimensionality == 2)>,
         typename = typename std::enable_if_t<(MultiArray2D::dimensionality==2)>,
         typename = void
        >
void Matrix2MAREF(char TA, MA const& A, MultiArray2D& M)
{
  using Type1 = typename std::decay<MA>::type::element;
  using Type2 = typename MultiArray2D::element;
  assert(TA=='N' || TA=='H' || TA=='T' || TA=='Z');
  if(TA=='N' || TA=='Z') {
    if(M.size(0) != A.size(0) or M.size(1) != A.size(1))
      throw std::runtime_error(" Error: Wrong dimensions in Matrix2MAREF.\n"); 
  } else if(TA=='T' || TA=='H') {
    if(M.size(0) != A.size(1) or M.size(1) != A.size(0))
      throw std::runtime_error(" Error: Wrong dimensions in Matrix2MAREF.\n");
  } else {
    throw std::runtime_error(" Error: Unknown operation in Matrix2MA.\n");
  }

  using ptrA = std::remove_cv_t<typename MA::element_ptr>;
  using ptrM = std::remove_cv_t<typename MultiArray2D::element_ptr>;

  if(TA == 'H') TA = 'C';
  if(TA=='Z') {
    // bad i gpu's
    for(int i=0; i<M.size(0); i++)
      for(int j=0; j<M.size(1); j++)
        M[i][j] = ma::conj(A[i][j]);
  } else if( not std::is_same<ptrA,ptrM>::value ) {
    if(TA=='N') {
      for(int i=0; i<M.size(0); i++)
        for(int j=0; j<M.size(1); j++)
          M[i][j] = A[i][j];
    } else if(TA=='T') {
      for(int i=0; i<M.size(0); i++)
        for(int j=0; j<M.size(1); j++)
          M[i][j] = A[j][i];
    } else if(TA=='C') {
      for(int i=0; i<M.size(0); i++)
        for(int j=0; j<M.size(1); j++)
          M[i][j] = ma::conj(A[j][i]);
    }
  } else {
    geam(   TA, TA, M.size(1), M.size(0),
          Type2(1.0), pointer_dispatch(A.origin()), A.stride(0),
          Type2(0.0), pointer_dispatch(A.origin()), A.stride(0),
          pointer_dispatch(M.origin()), M.stride(0)
    );
  }
}

template<class MA,
         class MultiArray2D,
         class Vector,
         typename = typename std::enable_if_t<(std::decay<MA>::type::dimensionality == 2)>,
         typename = typename std::enable_if_t<(MultiArray2D::dimensionality==2)>,
         typename = void
        >
void Matrix2MA(char TA, MA const& A, MultiArray2D& M, Vector const& occups)
{
  using Type1 = typename std::decay<MA>::type::element;
  using Type2 = typename MultiArray2D::element;
  assert(TA=='N' || TA=='H' || TA=='T' || TA=='Z');
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
  if(TA == 'H') TA = 'C';
  // bad i gpu's
  if(TA=='N') {
    for(int i=0; i<nrows; i++) {
      assert(occups[i] >= 0 && occups[i] < A.size(0));
      int ik = occups[i];
      for(int j=0; j<M.size(1); j++)
        M[i][j] = static_cast<Type2>(A[ik][j]);
    }
  } else if(TA=='T') {
    for(int i=0; i<nrows; i++) {
      assert(occups[i] >= 0 && occups[i] < A.size(0));
      int ik = occups[i];
      for(int j=0; j<M.size(1); j++)
        M[j][i] = static_cast<Type2>(A[ik][j]);
    }
  } else if(TA=='C') {
    for(int i=0; i<nrows; i++) {
      assert(occups[i] >= 0 && occups[i] < A.size(0));
      int ik = occups[i];
      for(int j=0; j<M.size(1); j++)
        M[j][i] = static_cast<Type2>(ma::conj(A[ik][j]));
    }
  } else if(TA=='Z') {
    for(int i=0; i<nrows; i++) {
      assert(occups[i] >= 0 && occups[i] < A.size(0));
      int ik = occups[i];
      for(int j=0; j<M.size(1); j++)
        M[i][j] = static_cast<Type2>(ma::conj(A[ik][j]));
    }
  }
}

}

#endif
