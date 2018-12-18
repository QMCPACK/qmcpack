////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
// Alfredo Correa, correaa@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////


#ifndef  AFQMC_APPLY_EXPM_HPP 
#define  AFQMC_APPLY_EXPM_HPP 

#include "AFQMC/Numerics/ma_operations.hpp"
#include <Utilities/FairDivide.h>

namespace qmcplusplus
{

namespace afqmc
{

namespace SlaterDeterminantOperations
{

namespace base
{

/*
 * Calculate S = exp(im*V)*S using a Taylor expansion of exp(V)
 */ 
template< class MatA,
          class MatB,
          class MatC
        >
inline void apply_expM( const MatA& V, MatB& S, MatC& T1, MatC& T2, int order=6)
{ 
  assert( V.shape()[0] == V.shape()[1] );
  assert( V.shape()[1] == S.shape()[0] );
  assert( S.shape()[0] == T1.shape()[0] );
  assert( S.shape()[1] == T1.shape()[1] );
  assert( S.shape()[0] == T2.shape()[0] );
  assert( S.shape()[1] == T2.shape()[1] );

  using ComplexType = typename std::decay<MatB>::type::element; 
  ComplexType zero(0.);
  MatC* pT1 = &T1;
  MatC* pT2 = &T2;

  T1 = S;
  for(int n=1; n<=order; n++) {
    ComplexType fact = ComplexType(0.0,1.0)*static_cast<ComplexType>(1.0/static_cast<double>(n));
    ma::product(fact,V,*pT1,zero,*pT2);
    // overload += ???
    for(int i=0, ie=S.shape()[0]; i<ie; i++)
     for(int j=0, je=S.shape()[1]; j<je; j++)
      S[i][j] += (*pT2)[i][j];
    std::swap(pT1,pT2);
  }

}

}

namespace shm 
{

/*
 * Calculate S = exp(im*V)*S using a Taylor expansion of exp(V)
 * V, S, T1, T2 are expected to be in shared memory.  
 */
template< class MatA,
          class MatB,
          class MatC,
          class communicator
        >
inline void apply_expM( const MatA& V, MatB& S, MatC& T1, MatC& T2, communicator& comm, int order=6)
{
  assert( V.shape()[0] == S.shape()[0] );
  assert( V.shape()[1] == S.shape()[0] );
  assert( S.shape()[0] == T1.shape()[0] );
  assert( S.shape()[1] == T1.shape()[1] );
  assert( S.shape()[0] == T2.shape()[0] );
  assert( S.shape()[1] == T2.shape()[1] );

  using boost::indices;
  using range_t = boost::multi_array_types::index_range;
  using ComplexType = typename std::decay<MatB>::type::element;

  const ComplexType zero(0.);
  const ComplexType im(0.0,1.0);
  MatC* pT1 = &T1;
  MatC* pT2 = &T2;

  int M0,Mn;
  std::tie(M0,Mn) = FairDivideBoundary(comm.rank(),int(S.shape()[0]),comm.size());

  assert( M0 <= Mn );  
  assert( M0 >= 0);

  T1[indices[range_t(M0,Mn)][range_t()]] = S[indices[range_t(M0,Mn)][range_t()]];
  comm.barrier();  
  for(int n=1; n<=order; n++) {
    const ComplexType fact = im*static_cast<ComplexType>(1.0/static_cast<double>(n));
    ma::product(fact,V[indices[range_t(M0,Mn)][range_t()]],*pT1,zero,(*pT2)[indices[range_t(M0,Mn)][range_t()]]);
    // overload += ???
    for(int i=M0; i<Mn; i++)
     for(int j=0, je=S.shape()[1]; j<je; j++)
      S[i][j] += (*pT2)[i][j];
    comm.barrier();  
    std::swap(pT1,pT2);
  }

}

}

}

}

}

#endif
