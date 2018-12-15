////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////


#ifndef  AFQMC_ROTATE_HPP 
#define  AFQMC_ROTATE_HPP 

#include <numeric>
#include "AFQMC/config.h"
#include <Utilities/FairDivide.h>
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/csr_blas.hpp"
#include "AFQMC/Matrix/matrix_emplace_wrapper.hpp"
//#include <boost/hana.hpp>
#include "AFQMC/Utilities/afqmc_TTI.hpp"
#include "mpi.h"

namespace qmcplusplus
{

namespace afqmc 
{

namespace sparse_rotate
{

/*
 *  Performs a (left) half rotation (and a possible transposition) of a Cholesky matrix. 
 *  The rotated matrix is sotored in "tuple" form in the provided container using emplace_back. 
 *  The generation of the matrix is distributed among the containers
 *  of all the cores in the node. 
 *  Input:
 *    -alpha: Sparse hermitian of rotation matrix for spin up, e.g. transpose(conj(Aup)).
 *    -beta: Sparse hermitian of rotation matrix for spin down, e.g. transpose(conj(Aup)).
 *    -A: Input sparse cholesky matrix.
 *    -cutoff: Value below which elements of rotated matrix are ignored.
 *  Output:
 *    -B: Container (e.g. std::vector) with a segment of the non-zero terms of the rotated matrix. 
 *        If transpose==false, the elements in the container are ordered by (a,k) values.
 *
 *  If transposed==true:
 *     B(n,ka-ka0) = sum_i^M alpha(a,i) * Spvn(ik,n) 
 *     B(n,ka+N*M-ka0) = sum_i^M beta(a,i) * Spvn(ki,n) 
 *  else:
 *     B(ka-ka0,n) = sum_i A(a,i) * Spvn(ik,n) 
 *     B(ka+N*M-ka0,n) = sum_i^M beta(a,i) * Spvn(ik,n), 
 *  where M/N is the number of rows/columns of alpha and beta.
 *  The number of rows of Spvn should be equal to M*M.
 *  If conjV=true, Spvn(ik,n) -> conj(Spvn(ki,n)) in the expressions above
 */ 
template< class Container,
          class task_group
        >
void halfRotateCholeskyMatrix(WALKER_TYPES type, task_group& TG, int k0, int kN, Container& Q, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, SpVType_shm_csr_matrix const& CholMat, bool transpose, bool conjV = false, double cutoff=1e-6, bool reserve_to_fit_=true)
{
  int NAEA = Alpha->shape()[0]; 
  int NAEB = Alpha->shape()[0]; 
  int NMO = Alpha->shape()[1]; 
  if(type==COLLINEAR)
    NAEB = Beta->shape()[0];
  int NEL = (type==CLOSED)?(NAEA):(NAEA+NAEB);
  int nvec = CholMat.shape()[1];  
  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();  

  assert(CholMat.shape()[0]==NMO*NMO); 
  assert(kN > k0); 
  if(type == CLOSED && kN > NMO)
    APP_ABORT(" Error: kN > NMO in halfRotateCholeskyMatrix. \n");

  // map from [0:2*NMO) to [0:NMO) in collinear case
  int k0_alpha, k0_beta, kN_alpha=0, kN_beta=0; 
  if(k0 >= 0) {
    k0_alpha = std::min(k0,NMO);
    kN_alpha = std::min(kN,NMO); 
    if(type==COLLINEAR) {
      kN_beta = std::max(kN,NMO)-NMO;
      k0_beta = std::max(k0,NMO)-NMO; 
    }
  } else {
    k0_alpha = k0_beta = kN_alpha = kN_beta = 0;
  }
   
  int ak0, ak1;
  int Qdim = NAEA*(kN_alpha-k0_alpha) + NAEB*(kN_beta-k0_beta);
  if(transpose) {
    if(not check_shape(Q,{nvec,Qdim}))
      APP_ABORT(" Error: Container Q has incorrect dimensions in halfRotateCholeskyMatrix. \n");
  } else {
    if(not check_shape(Q,{Qdim,nvec}))
      APP_ABORT(" Error: Container Q has incorrect dimensions in halfRotateCholeskyMatrix. \n");
  }  
  std::tie(ak0,ak1) = FairDivideBoundary(coreid,Qdim,ncores);

  if(type==NONCOLLINEAR)
    APP_ABORT(" GHF not yet implemented. \n");

  boost::multi_array<SPComplexType,1> vec(extents[nvec]);
  if(reserve_to_fit_) {
    std::vector<std::size_t> sz_per_row( Qdim ); 
    int cnt=0;
    for(int k=k0_alpha; k<kN_alpha; k++) {
      for(int a=0; a<NAEA; a++, cnt++) {
        if( cnt < ak0 ) continue;
        if( cnt >= ak1 ) break;  
        std::fill_n(vec.origin(),vec.size(),SPComplexType(0,0));
        auto Aa = (*Alpha)[a];
        for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i = Aa.non_zero_indices2_data()[ip];
          if(conjV)
            csr::axpy('C',Aai,CholMat[k*NMO+i],vec);  
          else
            csr::axpy('N',Aai,CholMat[i*NMO+k],vec);  
        }  
        if(transpose)  
          for(int n=0; n<nvec; n++)
            if(std::abs(vec[n])>cutoff) ++sz_per_row[n]; 
        else  
          for(int n=0; n<nvec; n++)
            if(std::abs(vec[n])>cutoff) ++sz_per_row[cnt];
      }
    }  
    // reset "amount of work done" to full alpha piece
    cnt=NAEA*(kN_alpha-k0_alpha);
    if(type==COLLINEAR) {
      // reset "shift"
      for(int k=k0_beta; k<kN_beta; k++) {
        for(int a=0; a<NAEB; a++, cnt++) {
          if( cnt < ak0 ) continue;
          if( cnt >= ak1 ) break;
          std::fill_n(vec.origin(),vec.size(),SPComplexType(0,0));
          auto Aa = (*Beta)[a];
          for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
            auto Aai = Aa.non_zero_values_data()[ip];
            auto i = Aa.non_zero_indices2_data()[ip];
            if(conjV)
              csr::axpy('C',Aai,CholMat[k*NMO+i],vec);
            else
              csr::axpy('N',Aai,CholMat[i*NMO+k],vec);
          }
          if(transpose)
            for(int n=0; n<nvec; n++)
              if(std::abs(vec[n])>cutoff) ++sz_per_row[n]; 
          else
            for(int n=0; n<nvec; n++)
              if(std::abs(vec[n])>cutoff) ++sz_per_row[cnt]; 
        }
      }
    }
    TG.Node().all_reduce_in_place_n(sz_per_row.begin(),sz_per_row.size(),std::plus<>());
    reserve_to_fit(Q,sz_per_row);
  }// else 
  //  Q.reserve( std::size_t(0.1 * nvec * NEL * NMO / ncores ) ); // by default, assume 10% sparsity
  
  int cnt=0;
  for(int k=k0_alpha; k<kN_alpha; k++) {
    for(int a=0; a<NAEA; a++, cnt++) {
      if( cnt < ak0 ) continue;
      if( cnt >= ak1 ) break;  
      std::fill_n(vec.origin(),vec.size(),SPComplexType(0,0));
      auto Aa = (*Alpha)[a];
      for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
        auto Aai = Aa.non_zero_values_data()[ip];
        auto i = Aa.non_zero_indices2_data()[ip];
        if(conjV)
          csr::axpy('C',Aai,CholMat[k*NMO+i],vec);
        else
          csr::axpy('N',Aai,CholMat[i*NMO+k],vec);
      }   
      if(transpose)  
        for(int n=0; n<nvec; n++)
          if(std::abs(vec[n])>cutoff) emplace(Q,std::forward_as_tuple(n,cnt,vec[n])); 
      else  
        for(int n=0; n<nvec; n++)
          if(std::abs(vec[n])>cutoff) emplace(Q,std::forward_as_tuple(cnt,n,vec[n])); 
    }
  }  
  // reset "amount of work done" to full alpha piece
  cnt=NAEA*(kN_alpha-k0_alpha);
  if(type==COLLINEAR) { 
    // reset "shift"
    for(int k=k0_beta; k<kN_beta; k++) {
      for(int a=0; a<NAEB; a++, cnt++) {
        if( cnt < ak0 ) continue;
        if( cnt >= ak1 ) break;
        std::fill_n(vec.origin(),vec.size(),SPComplexType(0,0));
        auto Aa = (*Beta)[a];
        for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i = Aa.non_zero_indices2_data()[ip];
          if(conjV)
            csr::axpy('C',Aai,CholMat[k*NMO+i],vec);
          else
            csr::axpy('N',Aai,CholMat[i*NMO+k],vec);
        }
        if(transpose)
          for(int n=0; n<nvec; n++)
            if(std::abs(vec[n])>cutoff) emplace(Q,std::forward_as_tuple(n,cnt,vec[n]));
        else
          for(int n=0; n<nvec; n++)
            if(std::abs(vec[n])>cutoff) emplace(Q,std::forward_as_tuple(cnt,n,vec[n]));
      }
    }
  }
}

/*
 * Calculates the rotated Cholesky matrix used in the calculation of the vias potential.
 *     v(n,ak) = sum_i A(a,i) * Spvn(ik,n)     A(a,i)=PsiT(i,a)*
 *     v(n,N*M+ak) = sum_i^M beta(a,i) * Spvn(ik,n), 
 *  where M/N is the number of rows/columns of alpha and beta.
 *  The number of rows of Spvn should be equal to M*M.
 *  Since we only rotate the "local" cholesky matrix, the algorithm is only parallelized
 *  over a node. In principle, it is possible to spread this over equivalent nodes.
 */ 
template<class task_group>
SpCType_shm_csr_matrix halfRotateCholeskyMatrixForBias(WALKER_TYPES type, task_group& TG, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, SpVType_shm_csr_matrix const& CholMat, double cutoff=1e-6)
{
  int NAEA = Alpha->shape()[0]; 
  int NAEB = Alpha->shape()[0]; 
  int NMO = Alpha->shape()[1]; 
  if(type!=CLOSED)
    NAEB = Beta->shape()[0];
  int NEL = (type==CLOSED)?(NAEA):(NAEA+NAEB);
  int nvec = CholMat.shape()[1];  
  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();  

// to speed up, generate new communicator for eqv_nodes and split full work among all
// cores in this comm. Then build from distributed container?

  assert(CholMat.shape()[0]==NMO*NMO); 

  std::size_t Qdim = NAEA*NMO;
  if(type==COLLINEAR) Qdim += NAEB*NMO;
  if(type==NONCOLLINEAR) Qdim = 2*NMO*(NAEA+NAEB);    
  std::size_t ak0, ak1;
  std::tie(ak0,ak1) = FairDivideBoundary(std::size_t(coreid),Qdim,std::size_t(ncores));

  if(type==NONCOLLINEAR)
    APP_ABORT(" GHF not yet implemented. \n");

  boost::multi_array<SPComplexType,1> vec(extents[nvec]);
  std::vector<std::size_t> sz_per_row( nvec ); 
  std::size_t cnt=0;
  for(int a=0; a<NAEA; a++) {
    for(int k=0; k<NMO; k++, cnt++) {
      if( cnt < ak0 ) continue;
      if( cnt >= ak1 ) break;  
      std::fill_n(vec.origin(),vec.size(),SPComplexType(0,0));
      auto Aa = (*Alpha)[a];
      for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
        auto Aai = Aa.non_zero_values_data()[ip];
        auto i = Aa.non_zero_indices2_data()[ip];
        csr::axpy('N',Aai,CholMat[i*NMO+k],vec);  
      }  
      for(std::size_t n=0; n<nvec; n++)
        if(std::abs(vec[n])>cutoff) ++(sz_per_row[n]); 
    }
  }  

  // reset "amount of work done" to full alpha piece
  cnt=NAEA*NMO;
  if(type==COLLINEAR) {
    // reset "shift"
    for(int a=0; a<NAEB; a++) {
      for(int k=0; k<NMO; k++, cnt++) {
        if( cnt < ak0 ) continue;
        if( cnt >= ak1 ) break;
        std::fill_n(vec.origin(),vec.size(),SPComplexType(0,0));
        auto Aa = (*Beta)[a];
        for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i = Aa.non_zero_indices2_data()[ip];
          csr::axpy('N',Aai,CholMat[i*NMO+k],vec);
        }
        for(std::size_t n=0; n<nvec; n++)
          if(std::abs(vec[n])>cutoff) ++sz_per_row[n]; 
      }
    }
  }
  TG.Node().all_reduce_in_place_n(sz_per_row.begin(),sz_per_row.size(),std::plus<>());

  using Alloc = boost::mpi3::intranode::allocator<SPComplexType>;
  SpCType_shm_csr_matrix::base ucsr({nvec,Qdim},{0,0},sz_per_row,Alloc(TG.Node()));

  using mat_wrapper = csr::matrix_emplace_wrapper<SpCType_shm_csr_matrix::base>;  
  mat_wrapper ucsr_wrapper(ucsr,TG.Node());

  cnt=0;
  for(int a=0; a<NAEA; a++) {
    for(int k=0; k<NMO; k++,cnt++) {
      if( cnt < ak0 ) continue;
      if( cnt >= ak1 ) break;  
      std::fill_n(vec.origin(),vec.size(),SPComplexType(0,0));
      auto Aa = (*Alpha)[a];
      for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
        auto Aai = Aa.non_zero_values_data()[ip];
        auto i = Aa.non_zero_indices2_data()[ip];
        csr::axpy('N',Aai,CholMat[i*NMO+k],vec);
      }   
      for(std::size_t n=0; n<nvec; n++)
        if(std::abs(vec[n])>cutoff) ucsr_wrapper.emplace(std::forward_as_tuple(n,cnt,vec[n])); 
    }
  }  
  // reset "amount of work done" to full alpha piece
  cnt=NAEA*NMO;
  if(type==COLLINEAR) { 
    // reset "shift"
    for(int a=0; a<NAEB; a++) {
      for(int k=0; k<NMO; k++, cnt++) {
        if( cnt < ak0 ) continue;
        if( cnt >= ak1 ) break;
        std::fill_n(vec.origin(),vec.size(),SPComplexType(0,0));
        auto Aa = (*Beta)[a];
        for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i = Aa.non_zero_indices2_data()[ip];
          csr::axpy('N',Aai,CholMat[i*NMO+k],vec);
        }
        for(std::size_t n=0; n<nvec; n++)
          if(std::abs(vec[n])>cutoff) ucsr_wrapper.emplace(std::forward_as_tuple(n,cnt,vec[n]));
      }
    }
  }
  ucsr_wrapper.push_buffer();
  return SpCType_shm_csr_matrix(std::move(ucsr)); 
}

}

namespace ma_rotate
{

template< class MultiArray2D, 
          class task_group
        >
void halfRotateCholeskyMatrix(WALKER_TYPES type, task_group& TG, int k0, int kN, MultiArray2D&& Q, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, SpVType_shm_csr_matrix const& CholMat, bool transpose, bool conjV = false, double cutoff=1e-6)
{
  int NAEA = Alpha->shape()[0]; 
  int NAEB = 0; 
  int NMO = Alpha->shape()[1]; 
  if(type==COLLINEAR)
    NAEB = Beta->shape()[0];
  int NEL = (type==CLOSED)?(NAEA):(NAEA+NAEB);
  int nvec = CholMat.shape()[1];  
  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();  

  assert(CholMat.shape()[0]==NMO*NMO); 
  if(type == CLOSED && kN > NMO)
    APP_ABORT(" Error: kN > NMO in halfRotateCholeskyMatrix. \n");

  // map from [0:2*NMO) to [0:NMO) in collinear case
  int k0_alpha, k0_beta, kN_alpha=0, kN_beta=0;
  if(k0 >= 0) {
    k0_alpha = std::min(k0,NMO);
    kN_alpha = std::min(kN,NMO);
    if(type==COLLINEAR) {
      k0_beta = std::max(k0,NMO)-NMO;
      kN_beta = std::max(kN,NMO)-NMO;
    }
  } else {
    k0_alpha = k0_beta = kN_alpha = kN_beta = 0;
  }

  int ak0, ak1;
  int Qdim = NAEA*(kN_alpha-k0_alpha) + NAEB*(kN_beta-k0_beta);
  if(transpose) {
    assert(Q.shape()[0]==nvec);
    assert(Q.shape()[1]==Qdim);
  } else {
    assert(Q.shape()[0]==Qdim);
    assert(Q.shape()[1]==nvec);
  }
  std::tie(ak0,ak1) = FairDivideBoundary(coreid,Qdim,ncores);

  if(type==NONCOLLINEAR)
    APP_ABORT(" GHF not yet implemented. \n");

  int cnt=0;
  for(int k=k0_alpha; k<kN_alpha; k++) {
    for(int a=0; a<NAEA; a++, cnt++) {
      if( cnt < ak0 ) continue;
      if( cnt >= ak1 ) break;  
      if(transpose) {
        auto vec = Q[indices[range_t()][cnt]];
        for(auto& v:vec) v=SPComplexType(0,0); 
        auto Aa = (*Alpha)[a];
        for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i = Aa.non_zero_indices2_data()[ip];
          if(conjV)
            csr::axpy('C',Aai,CholMat[k*NMO+i],vec);
          else
            csr::axpy('N',Aai,CholMat[i*NMO+k],vec);
        }   
      } else {
        auto vec = Q[cnt];
        for(auto& v:vec) v=SPComplexType(0,0);
        auto Aa = (*Alpha)[a];
        for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i = Aa.non_zero_indices2_data()[ip];
          if(conjV)
            csr::axpy('C',Aai,CholMat[k*NMO+i],vec);
          else
            csr::axpy('N',Aai,CholMat[i*NMO+k],vec);
        }
      }
    }
  }  
  // reset "amount of work done" to full alpha piece
  cnt=NAEA*(kN_alpha-k0_alpha);
  if(type==COLLINEAR) {
    // reset "shift"
    for(int k=k0_beta; k<kN_beta; k++) {
      for(int a=0; a<NAEB; a++, cnt++) {
        if( cnt < ak0 ) continue;
        if( cnt >= ak1 ) break;
        if(transpose) {
          auto vec = Q[indices[range_t()][cnt]];
          for(auto& v:vec) v=SPComplexType(0,0);
          auto Aa = (*Beta)[a];
          for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
            auto Aai = Aa.non_zero_values_data()[ip];
            auto i = Aa.non_zero_indices2_data()[ip];
            if(conjV)
              csr::axpy('C',Aai,CholMat[k*NMO+i],vec);
            else
              csr::axpy('N',Aai,CholMat[i*NMO+k],vec);
          }
        } else {
          auto vec = Q[cnt];
          for(auto& v:vec) v=SPComplexType(0,0);
          auto Aa = (*Beta)[a];
          for(int ip = 0; ip<Aa.num_non_zero_elements(); ++ip) {
            auto Aai = Aa.non_zero_values_data()[ip];
            auto i = Aa.non_zero_indices2_data()[ip];
            if(conjV)
              csr::axpy('C',Aai,CholMat[k*NMO+i],vec);
            else
              csr::axpy('N',Aai,CholMat[i*NMO+k],vec);
          }
        }

      }
    }
  }
  TG.Node().barrier();
}

}

}

}

#endif
