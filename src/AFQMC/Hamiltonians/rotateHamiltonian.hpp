//////////////////////////////////////////////////////////////////////////////
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

#ifndef QMCPLUSPLUS_AFQMC_ROTATEHAMILTONIAN_HPP
#define QMCPLUSPLUS_AFQMC_ROTATEHAMILTONIAN_HPP

#include<vector>
#include<tuple>
#include<mpi.h>
#include<algorithm>
#include<numeric>

#include "Numerics/OhmmsBlas.h"
#include "AFQMC/config.h"
#include "AFQMC/Numerics/DenseMatrixOperations.h"
#include "AFQMC/Numerics/SparseMatrixOperations.h"
#include "AFQMC/Utilities/taskgroup.h"
//#include "AFQMC/Matrix/spma_communications.hpp"
#include "AFQMC/Matrix/mpi3_SHMBuffer.hpp"
#include "AFQMC/Matrix/matrix_emplace_wrapper.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/SlaterDeterminantOperations/rotate.hpp"
#include "AFQMC/Utilities/afqmc_TTI.hpp"

#include "AFQMC/Hamiltonians/rotateHamiltonian_Helper2.hpp"
#include "AFQMC/Hamiltonians/rotateHamiltonian_Helper3.hpp"

namespace qmcplusplus
{

namespace afqmc
{

inline void check_wavefunction_consistency(WALKER_TYPES type, PsiT_Matrix *A, PsiT_Matrix *B, int NMO, int NAEA, int NAEB) 
{
    if(type == CLOSED) {
      if(A->shape()[1] != NMO || A->shape()[0] != NAEA) {
        app_error()<<" Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency(): wfn_type=0, NMO, NAEA, A.rows, A.cols: " <<NMO <<" " <<NAEA <<" " <<A->shape()[0] <<" " <<A->shape()[1] <<std::endl; 
        APP_ABORT(" Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency().\n");
      }
    } else if(type == COLLINEAR) {
      if(A->shape()[1] != NMO || A->shape()[0] != NAEA || B->shape()[1] != NMO || B->shape()[0] != NAEB) {
        app_error()<<" Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency(): wfn_type=1, NMO, NAEA, NAEB, A.rows, A.cols, B.rows, B.cols: " 
        <<NMO <<" " <<NAEA <<" " <<NAEB <<" " 
        <<A->shape()[0] <<" " <<A->shape()[1] <<" "
        <<B->shape()[0] <<" " <<B->shape()[1] <<std::endl; 
        APP_ABORT(" Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency().\n");
      }
    } else if(type==NONCOLLINEAR) {
      if(A->shape()[1] != 2*NMO || A->shape()[0] != (NAEB+NAEA)) {
        app_error()<<" Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency(): wfn_type=1, NMO, NAEA, NAEB, A.rows, A.cols: " <<NMO <<" " <<NAEA <<" " <<NAEB <<" " <<A->shape()[0] <<" " <<A->shape()[1] <<std::endl; 
        APP_ABORT(" Error: Incorrect Slater Matrix dimensions in check_wavefunction_consistency().\n");
      }
    } else {
      app_error()<<" Error: Unacceptable walker_type in check_wavefunction_consistency(): " <<type <<std::endl;
      APP_ABORT(" Error: Unacceptable walker_type in check_wavefunction_consistency(). \n");
    }
}

inline boost::multi_array<SPComplexType,1> rotateHij(WALKER_TYPES walker_type, int NMO, int NAEA, int NAEB, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, const std::vector<s2D<ValueType> >& H1)
{
  boost::multi_array<SPComplexType,1> N;  
  boost::multi_array<ComplexType,2> M;  
  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0);

  // 1-body part 
  if(walker_type == CLOSED) {

    ValueType V;
    M.resize(extents[NMO][NMO]);
    std::fill_n(M.origin(),NMO*NMO,ComplexType(0.0,0.0));
    N.resize(extents[NAEA*NMO]);
#if(AFQMC_SP)
    boost::multi_array<ComplexType,2> N_(extents[NAEA][NMO]);
#else
    boost::multi_array_ref<ComplexType,2> N_(N.origin(),extents[NAEA][NMO]);
#endif

    std::vector<s2D<ValueType> >::const_iterator it1 = H1.begin();
    while(it1 != H1.end()) {
      IndexType i,j;
      std::tie (i,j,V) = *it1++;
      M[i][j] = ValueType(2.0)*V; //
      if( i!=j ) M[j][i] = ValueType(2.0)*myconj(V);
    }
    
    ma::product(*Alpha,M,N_);
#if(AFQMC_SP)
    std::copy_n(N_.origin(),NAEA*NMO,N.origin());
#endif

  } else if(walker_type == COLLINEAR) {

    ValueType V;
    M.resize(extents[NMO][NMO]);
    std::fill_n(M.origin(),NMO*NMO,ComplexType(0.0,0.0));
    N.resize(extents[(NAEA+NAEB)*NMO]);
#if(AFQMC_SP)
    boost::multi_array<ComplexType,2> NA_(extents[NAEA][NMO]);
    boost::multi_array<ComplexType,2> NB_(extents[NAEB][NMO]);
#else
    boost::multi_array_ref<ComplexType,2> NA_(N.origin(),extents[NAEA][NMO]);
    boost::multi_array_ref<ComplexType,2> NB_(N.origin()+NAEA*NMO,extents[NAEB][NMO]);
#endif

    std::vector<s2D<ValueType> >::const_iterator it1 = H1.begin();
    while(it1 != H1.end()) {
      IndexType i,j;
      std::tie (i,j,V) = *it1++;
      if(i<NMO && j<NMO) {
        M[i][j] = V; //
        if( i!=j ) M[j][i] = myconj(V);
      }
    }

    ma::product(*Alpha,M,NA_);
    ma::product(*Beta,M,NB_);
#if(AFQMC_SP)
    std::copy_n(NA_.origin(),NAEA*NMO,N.origin());
    std::copy_n(NB_.origin(),NAEB*NMO,N.origin()+NAEA*NMO);
#endif

  } else if(walker_type == NONCOLLINEAR) {

    ValueType V;
    M.resize(extents[2*NMO][2*NMO]);
    std::fill_n(M.origin(),4*NMO*NMO,ComplexType(0.0,0.0));
    N.resize(extents[(NAEA+NAEB)*2*NMO]);
#if(AFQMC_SP)
    boost::multi_array<ComplexType,2> N_(extents[NAEA+NAEB][2*NMO]);
#else
    boost::multi_array_ref<ComplexType,2> N_(N.origin(),extents[NAEA+NAEB][2*NMO]);
#endif

    std::vector<s2D<ValueType> >::const_iterator it1 = H1.begin();
    while(it1 != H1.end()) {
      IndexType i,j;
      std::tie (i,j,V) = *it1++;
      M[i][j] = V; //
      if( i!=j ) M[j][i] = myconj(V);
    }

    ma::product(*Alpha,M,N_);
#if(AFQMC_SP)
    std::copy_n(N_.origin(),(NAEA+NAEB)*2*NMO,N.origin());
#endif

  }

  return N;
}

inline boost::multi_array<SPComplexType,1> rotateHij(WALKER_TYPES walker_type, int NMO, int NAEA, int NAEB, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, const boost::multi_array<ComplexType,2>& H1)
{
  boost::multi_array<SPComplexType,1> N;  
  const ComplexType one = ComplexType(1.0);
  const ComplexType zero = ComplexType(0.0);

  // 1-body part 
  if(walker_type == CLOSED) {

    N.resize(extents[NAEA*NMO]);
#if(AFQMC_SP)
    boost::multi_array<ComplexType,2> N_(extents[NAEA][NMO]);
#else
    boost::multi_array_ref<ComplexType,2> N_(N.origin(),extents[NAEA][NMO]);
#endif

    ma::product(*Alpha,H1,N_);
#if(AFQMC_SP)
    std::copy_n(N_.origin(),NAEA*NMO,N.origin());
#endif
    ma::scal(SPComplexType(2.0),N);  

  } else if(walker_type == COLLINEAR) {

    N.resize(extents[(NAEA+NAEB)*NMO]);
#if(AFQMC_SP)
    boost::multi_array<ComplexType,2> NA_(extents[NAEA][NMO]);
    boost::multi_array<ComplexType,2> NB_(extents[NAEB][NMO]);
#else
    boost::multi_array_ref<ComplexType,2> NA_(N.origin(),extents[NAEA][NMO]);
    boost::multi_array_ref<ComplexType,2> NB_(N.origin()+NAEA*NMO,extents[NAEB][NMO]);
#endif

    ma::product(*Alpha,H1,NA_);
    ma::product(*Beta,H1,NB_);
#if(AFQMC_SP)
    std::copy_n(NA_.origin(),NAEA*NMO,N.origin());
    std::copy_n(NB_.origin(),NAEB*NMO,N.origin()+NAEA*NMO);
#endif

  } else if(walker_type == NONCOLLINEAR) {

    N.resize(extents[(NAEA+NAEB)*2*NMO]);
#if(AFQMC_SP)
    boost::multi_array<ComplexType,2> N_(extents[NAEA+NAEB][2*NMO]);
#else
    boost::multi_array_ref<ComplexType,2> N_(N.origin(),extents[NAEA+NAEB][2*NMO]);
#endif

    ma::product(*Alpha,H1,N_);
#if(AFQMC_SP)
    std::copy_n(N_.origin(),(NAEA+NAEB)*2*NMO,N.origin());
#endif

  }

  return N;
}

template<class Container = std::vector<std::tuple<int,int,SPComplexType>>>
inline void rotateHijkl(std::string& type, WALKER_TYPES walker_type, TaskGroup_& TG, Container& Vijkl, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, SpVType_shm_csr_matrix const& V2_fact, const RealType cut, int maximum_buffer_size, bool reserve_to_fit_=true, bool global_reserve = true)
{
  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  bool distribute_Ham  = TG.getNumberOfTGs() > 1;
  if(distribute_Ham)
    APP_ABORT(" Distributed V2_fact not yet implemented. \n");

  int NAEA = Alpha->shape()[0];
  int NMO = Alpha->shape()[1];
  int NAEB = NAEA;
  if( walker_type == COLLINEAR ) NAEB = Beta->shape()[0];

  // <ab||kl> = sum_n Qk(a,n) * Rl(b,n) - Rl(a,n)*Qk(b,n),
  // where:
  //   Qk(a,n) = sum_i conj(Amat(i,a)) * V2_fact(ik,n)
  //   Rl(a,n) = sum_i conj(Amat(i,a)) * conj(V2_fact(li,n))
  // For real build, Qk=Rk
  //
  // For parallelization, distribute (k,l) pairs over nodes.
  // Build ahead of time Qk/Rl matrices in shared memory to reduce memory/setup time.
  // Assemble integrals in parallel and fully distributed.
  // Collect on all nodes.
  //    - For distributed hamiltonians, you do not need to worry about keeping contiguous
  //    segments of the hamiltonian. Only that the distribution over nodes is roughly equal.
  //    Write a simple algorithm that balances the number of terms in a TG. 
  //  

  bool sparseQk = (type == "SD" || type == "SS");
  bool sparseRl = (type == "SS");

  app_log()<<" Calculating half-rotated Hamiltonian using ";
  if(type == "SD")
    app_log()<<"Sparse x Dense";
  else if(type == "SS")
    app_log()<<"Sparse x Sparse";
  else if(type == "DD")
    app_log()<<"Dense x Dense";
  app_log()<<" matrix multiplication. \n";

  mpi3_SHMBuffer<SPComplexType> tQk_shmbuff(TG.Node(),1);
  mpi3_SHMBuffer<SPComplexType> Qk_shmbuff(TG.Node(),1);
  mpi3_SHMBuffer<SPComplexType> Rl_shmbuff(TG.Node(),1);

  int NMO2 = (walker_type == CLOSED)?NMO:2*NMO;
  int ngrp = std::min(NMO2,nnodes);
  std::vector<int> M_split(ngrp+1);

  // split the orbitals among processors:  ngrp/2 + ngrp%2 for alpha, ngrp/2 for beta 
  M_split[0]=0;
  if(walker_type == CLOSED) {
    FairDivide(NMO2,ngrp,M_split);
  } else if (walker_type == COLLINEAR) {
    if(ngrp==1)
      APP_ABORT(" Error: Current implementation of rotateHijkl requires at least 2 nodes.\n");
    // This should be 2/3-1/3 partitioning between alpha/beta
    // to balance the resulting matrix elements better
    int nalpha = ngrp/2;
    for(int i=0; i<nalpha; i++) {
      int m0,m1;
      std::tie(m0,m1) = FairDivideBoundary(i,NMO,nalpha);
      M_split[i+1]=m1;
    }
    assert(M_split[nalpha]==NMO);
    for(int i=0, nbeta=ngrp/2+ngrp%2; i<nbeta; i++) {
      int m0,m1;
      std::tie(m0,m1) = FairDivideBoundary(i,NMO,nbeta);
      M_split[i+nalpha+1]=NMO+m1;
    }
    assert(M_split[ngrp]==2*NMO);
  } else if (walker_type == NONCOLLINEAR) {
    APP_ABORT(" Finish. \n");
  }

  // Construct your set of Q(k,a,m), R(l,b,m)
  int l0 = (nodeid<ngrp)?(M_split[nodeid]):(-1);
  int lN = (nodeid<ngrp)?(M_split[nodeid+1]):(-1);
  bool amIAlpha = true; // for simplicity, subset of bands must be either elpha or beta, not mixed
  if( l0 < NMO && (lN-1) < NMO )
    amIAlpha = true;
  else if( l0 >= NMO && lN >= NMO )
    amIAlpha = false;
  else {
    std::cerr<<"l0, lN, nodeid, ngrp, NMO: " <<l0 <<" " <<lN <<" " <<nodeid <<" " <<ngrp <<" " <<NMO <<std::endl;
    std::cerr<<" Error: Current algorithm requires an even number of processors. \n\n\n";
    APP_ABORT(" Error: Current algorithm requires an even number of processors. \n\n\n");
  }
  int norb = lN-l0;
  int maxnorb = 0;
  for(int i=0; i<ngrp; i++) maxnorb = std::max(maxnorb,M_split[i+1]-M_split[i]);
  int nvec = V2_fact.shape()[1];
  // must gather over heads of TG to get nchol per TG and total # chol vecs 
  // Rl(k, a, m), k:[0:NMO2], a:[0:NAEA], m:[0:nvec] 

  app_log()<<" Approximate memory usage for half-rotated Hamiltonian construction: \n"
           <<"   max. number of orbital in a node: " <<maxnorb <<"\n"
           <<"   Qk/Rl matrices size (assuming dense) each = maxnorb * nup * ncholvecs complex numbers = "
                    <<sizeof(SPComplexType)*maxnorb*NAEA*nvec/1024.0/1024.0 <<" MB \n"
           <<"   Maximum size of communication buffer: " <<maximum_buffer_size <<" MB" <<std::endl;

  const int nrow = norb * ((amIAlpha)?NAEA:NAEB);
  const int ncol = nvec;
  int dummy_nrow=nrow, dummy_ncol=ncol;
  int mat_size = nrow*ncol; 
  if(nodeid < ngrp) {
    if(not sparseQk) {
      Qk_shmbuff.resize(mat_size);
      if(coreid==0) std::fill_n(Qk_shmbuff.data(),Qk_shmbuff.size(),SPComplexType(0.0));
    } 
    if(not sparseRl) {
      Rl_shmbuff.resize(mat_size);
      if(coreid==0) std::fill_n(Rl_shmbuff.data(),Rl_shmbuff.size(),SPComplexType(0.0));
    } 
  } else {
    dummy_nrow=dummy_ncol=0;
  }

  using Alloc = boost::mpi3::intranode::allocator<SPComplexType>;
  // global_origin is not set correctly, careful not to rely on it
  SpCType_shm_csr_matrix SpQk({nrow,ncol},{0,0},0,Alloc(TG.Node()));
  SpCType_shm_csr_matrix SpRl({ncol,nrow},{0,0},0,Alloc(TG.Node()));

  if(sparseQk)  dummy_nrow=dummy_ncol=0;
  boost::multi_array_ref<SPComplexType,2> Qk(Qk_shmbuff.data(),extents[dummy_nrow][dummy_ncol]);  
  dummy_nrow=nrow; dummy_ncol=ncol;
  if(sparseRl or nodeid >= ngrp )  dummy_nrow=dummy_ncol=0;
  boost::multi_array_ref<SPComplexType,2> Rl(Rl_shmbuff.data(),extents[dummy_ncol][dummy_nrow]);  

  if(distribute_Ham) {
   APP_ABORT(" Finish THIS (43)!!! \n\n\n");
  } else {

    //   Q(k,a,n) = sum_i conj(Amat(i,a)) * V2_fact(ik,n)
    //   R(l,a,n) = sum_i conj(Amat(i,a)) * conj(V2_fact(li,n))

    // Construct SpQk[k,n,nvec]
    if(sparseQk) {
      sparse_rotate::halfRotateCholeskyMatrix(walker_type,TG,l0,lN,SpQk,Alpha,Beta,V2_fact,false,false,cut,true);
      SpQk.remove_empty_spaces();  // just in case  
    } else
      ma_rotate::halfRotateCholeskyMatrix(walker_type,TG,l0,lN,Qk,Alpha,Beta,V2_fact,false,false,cut);

#if defined(QMC_COMPLEX)
   // Construct SpRl[nvec,k,a] 
    if(sparseRl) {
      // since Rl is transposed, I can't emplace_back on the csr matrix. 
      // In this case, I need to use a temporary ucsr with an emplace_wrapper
      using ucsr_matrix = ma::sparse::ucsr_matrix<SPComplexType,int,std::size_t,
                                boost::mpi3::intranode::allocator<SPComplexType>,
                                boost::mpi3::intranode::is_root>;
      ucsr_matrix ucsr({ncol,nrow},{0,0},0,Alloc(TG.Node()));
      csr::matrix_emplace_wrapper<ucsr_matrix> ucsr_wrapper(ucsr,TG.Node()); 
      sparse_rotate::halfRotateCholeskyMatrix(walker_type,TG,l0,lN,ucsr_wrapper,Alpha,Beta,V2_fact,true,true,cut,true);  
      ucsr_wrapper.push_buffer(); // push any remaining elements in temporary buffer
      SpRl = std::move(ucsr);
      SpRl.remove_empty_spaces();  // just in case  
    } else
      ma_rotate::halfRotateCholeskyMatrix(walker_type,TG,l0,lN,Rl,Alpha,Beta,V2_fact,true,true,cut);
#else
    if(sparseRl) {
      if(sparseQk) {
        SpRl = std::move(csr::shm::transpose(SpQk));
        SpRl.remove_empty_spaces();  // just in case  
      } else {
        app_error()<<" Error: Incorrect matrix setup in createHamiltonianForGeneralDeterminant. sparseRl=True, sparseQk=False."
                   <<std::endl;
        APP_ABORT(" Error: Incorrect matrix setup in createHamiltonianForGeneralDeterminant. sparseRl=True, sparseQk=False. \n");
      }
    } else {
      if(sparseQk) {
        csr::shm::transpose(SpQk,Rl);  
      } else {
        // Qk[norb*NAEA,nvec]
        // Rl[nvec,norb*NAEA]
        int n0_,n1_,sz_ = Qk.shape()[0];
        std::tie(n0_, n1_) = FairDivideBoundary(coreid,sz_,ncores);
        if(n1_-n0_>0) 
          ma::tranpose(Qk[indices[range_t(n0_,n1_)][range_t()]],Rl[range_t()][range_t(n0_,n1_)]); 
      }
    }
#endif
  }

  TG.node_barrier();
  // let the maximum message be maximum_buffer_size MB
  int maxnt = std::max(1,static_cast<int>(std::floor(maximum_buffer_size*1024.0*1024.0/sizeof(SPComplexType))));

  // control the size of the MPI messages.
  // communicate blocks of Qk up to maxnt terms  
  std::vector<int> nkbounds;          // local bounds for communication  
  std::vector<int> Qknum(nnodes);     // number of blocks per node
  std::vector<int> Qksizes;           // number of terms and number of k vectors in block for all nodes
  if(coreid==0) {
    int n_=0, ntcnt=0, n0=0;
    int NEL = (amIAlpha)?NAEA:NAEB;
    for(int i=0; i<norb; i++) {
      int ntt;
      if(sparseQk)
        ntt = SpQk.num_non_zero_elements(i); 
      else
        ntt = nvec*NAEA;
      assert(ntt < maxnt);
      if(ntcnt+ntt > maxnt) {
        nkbounds.push_back(ntcnt);
        nkbounds.push_back(i-n0);
        ntcnt=ntt;
        n0=i;
        n_++;
      } else {
        ntcnt+=ntt;
      }
    }
    if(ntcnt > 0) {
      // push last block
      n_++;
      nkbounds.push_back(ntcnt);
      nkbounds.push_back(norb-n0);
    }

    MPI_Allgather(&n_,1,MPI_INT,Qknum.data(),1,MPI_INT,TG.Cores().impl_);

    int ntt = std::accumulate(Qknum.begin(),Qknum.end(),0);
    Qksizes.resize(2*ntt);

    std::vector<int> cnts(nnodes);
    std::vector<int> disp(nnodes);
    int cnt=0;
    for(int i=0; i<nnodes; i++) {
      cnts[i] = Qknum[i]*2;
      disp[i]=cnt;
      cnt+=cnts[i];
    }
    MPI_Allgatherv(nkbounds.data(),nkbounds.size(),MPI_INT,Qksizes.data(),cnts.data(),disp.data(),MPI_INT,TG.Cores().impl_ );

  }

  MPI_Bcast(Qknum.data(),nnodes,MPI_INT,0,TG.Node().impl_);
  int ntt = std::accumulate(Qknum.begin(),Qknum.end(),0);
  if(!coreid==0)
    Qksizes.resize(2*ntt);
  MPI_Bcast(Qksizes.data(),Qksizes.size(),MPI_INT,0,TG.Node().impl_);

// store {nterms,nk} for all nodes 
// use it to know communication pattern 

  int maxnk = 0;        // maximum number of k vectors in communication block 
  long maxqksize = 0;   // maximum size of communication block
  for(int i=0; i<ntt; i++) {
    if(Qksizes[2*i] > maxqksize) maxqksize = Qksizes[2*i];
    if(Qksizes[2*i+1] > maxnk) maxnk = Qksizes[2*i+1];
  }

  app_log()<<"   Maximum number of (k,l) pairs communicated in a block: " <<maxnk <<"\n"
           <<"   Temporary integral matrix Ta: " <<norb*NAEA*maxnk*NAEA*sizeof(SPComplexType)/1024.0/1024.0 <<" MB " <<std::endl;

  // temporary shared memory space for local "dense" result
  mpi3_SHMBuffer<SPComplexType> Ta_shmbuff(TG.Node(),norb*NAEA*maxnk*NAEA);

  // setup working sparse matrix  
  dummy_nrow=maxnk * NAEA; dummy_ncol=nvec;
  SpCType_shm_csr_matrix SptQk({maxnk * NAEA,nvec},{0,0},0,Alloc(TG.Node()));
  if(sparseQk) { 
    std::size_t sz_ = std::ceil(maxqksize/SptQk.shape()[0]);
    SptQk.reserve(sz_);
  } else 
    tQk_shmbuff.resize(maxnk * NAEA * nvec);
  if(sparseQk)  dummy_nrow=dummy_ncol=0;

  myTimer Timer_;

  if(reserve_to_fit_) {
    // count and resize container
    std::vector<std::size_t> sz_local;
    if(walker_type==CLOSED) sz_local.resize(NMO*NAEA);
    else if(walker_type==COLLINEAR) sz_local.resize(NMO*(NAEA+NAEB));
    else if(walker_type==NONCOLLINEAR) sz_local.resize(2*NMO*(NAEA+NAEB));
    
    for(int nn=0, nb=0, nkcum=0; nn<ngrp; nn++) {

      // just checking
      assert(nkcum==M_split[nn]);
      if(M_split[nn+1]==M_split[nn]) continue;
      int nblk = Qknum[nn];
      long ntermscum=0;
      for( int bi = 0; bi < nblk; bi++, nb++) {
        int nterms = Qksizes[2*nb];      // number of terms in block 
        int nk = Qksizes[2*nb+1];        // number of k-blocks in block
        int k0 = nkcum;                  // first value of k in block
        nkcum+=nk;
        int kN = nkcum;                  // last+1 value
        int NEL0 = (k0<NMO)?NAEA:NAEB;   // number of electrons in this spin block
        assert(nk > 0 && nk <= maxnk );  // just checking

        boost::multi_array_ref<SPComplexType,2> tQk(tQk_shmbuff.data(),extents[nk*NEL0][nvec]);

        Timer_.reset("T0");
        Timer_.start("T0");
        if(sparseQk) {
          if(coreid==0) {
            if(nn == nodeid) {
              auto ka0 = (k0-M_split[nn])*NEL0;
              auto kaN = (k0-M_split[nn]+nk)*NEL0;
              auto n0 = *SpQk.pointers_begin( ka0 );
              auto n1 = *SpQk.pointers_end(kaN);
              int nt_ = static_cast<int>(n1-n0);
              assert(ntermscum==n0);
              assert(nt_==nterms);
              std::copy(std::addressof(*SpQk.non_zero_values_data(ka0)),
                        std::addressof(*SpQk.non_zero_values_data(kaN)),
                        std::addressof(*SptQk.non_zero_values_data()));  
              std::copy(std::addressof(*SpQk.non_zero_indices2_data(ka0)),
                        std::addressof(*SpQk.non_zero_indices2_data(kaN)),
                        std::addressof(*SptQk.non_zero_indices2_data()));  
              for(int i=0, j=ka0; i<=nk*NEL0; i++, j++) {
                SptQk.pointers_begin()[i] = SpQk.pointers_begin()[j]-n0;
                SptQk.pointers_end()[i] = SpQk.pointers_end()[j]-n0;
              }
            }
            TG.Cores().broadcast_value(nterms,nn);  
            TG.Cores().broadcast_n(std::addressof(*SptQk.pointers_begin()),
                                  SptQk.shape()[0],nn);
            TG.Cores().broadcast_n(std::addressof(*SptQk.pointers_end()),
                                  SptQk.shape()[0],nn);
            TG.Cores().broadcast_n(std::addressof(*SptQk.non_zero_values_data()),
                                  nterms,nn);  
            TG.Cores().broadcast_n(std::addressof(*SptQk.non_zero_indices2_data()),
                                  nterms,nn);  
          }
          TG.node_barrier();
          // for safety, keep track of sum
          ntermscum += static_cast<long>(nterms);
        } else {
          if(coreid==0) {
            if(nn == nodeid)
              std::copy(Qk.origin()+bi*maxnk*NEL0*nvec,Qk.origin()+(bi*maxnk+nk)*NEL0*nvec,tQk.origin());
            TG.Cores().broadcast_n(tQk.origin(),nk*NEL0*nvec,nn);
          }
          TG.node_barrier();
        }
        Timer_.stop("T0");
        app_log()<<" Loop: " <<nn <<"/" <<ngrp <<" " <<bi <<"/" <<nblk 
                 <<" communication: " <<Timer_.total("T0") <<" "; 

        boost::multi_array_ref<ComplexType,2> Ta(Ta_shmbuff.data(),extents[nk*NEL0][nrow]);

        Timer_.reset("T0");
        Timer_.start("T0");
        if(type == "SD")
          count_Qk_x_Rl(walker_type,TG,sz_local,k0,kN,l0,lN,NMO,NAEA,NAEB,SptQk[{0,std::size_t(nk*NEL0)}],Rl,Ta,cut);
        else if(type == "DD")   
          count_Qk_x_Rl(walker_type,TG,sz_local,k0,kN,l0,lN,NMO,NAEA,NAEB,tQk,Rl,Ta,cut);
        Timer_.stop("T0");
        app_log()<<" QxR: " <<Timer_.total("T0") <<std::endl; 
//      else if(type == "SS")

      }
    }

    std::size_t tot_sz_local = std::accumulate(sz_local.begin(),sz_local.end(),std::size_t(0));

    std::vector<std::size_t> sz_global(sz_local.size());
    TG.Global().all_reduce_n(sz_local.begin(),sz_local.size(),sz_global.begin(),std::plus<>());
    std::size_t tot_sz_global = std::accumulate(sz_global.begin(),sz_global.end(),std::size_t(0));

    std::vector<std::size_t> sz_node(sz_local.size());
    TG.Node().all_reduce_n(sz_local.begin(),sz_local.size(),sz_node.begin(),std::plus<>());
    std::size_t tot_sz_node = std::accumulate(sz_node.begin(),sz_node.end(),std::size_t(0));

    std::size_t sz_node_min= TG.Global().all_reduce_value(tot_sz_node,boost::mpi3::min<>());
    std::size_t sz_node_max= TG.Global().all_reduce_value(tot_sz_node,boost::mpi3::max<>());
    std::size_t sz_local_min= TG.Global().all_reduce_value(tot_sz_local,boost::mpi3::min<>());
    std::size_t sz_local_max= TG.Global().all_reduce_value(tot_sz_local,boost::mpi3::max<>());

    app_log()<<"  Number of terms in Vijkl: \n" 
           <<"    Local: (min/max)" 
           <<sz_local_min <<" " 
           <<sz_local_min*(sizeof(SPComplexType)+sizeof(int))/1024.0/1024.0 <<" MB  -  " 
           <<sz_local_max <<" " 
           <<sz_local_max*(sizeof(SPComplexType)+sizeof(int))/1024.0/1024.0 <<" MB \n" 
           <<"    Node (min/max): " 
           <<sz_node_min <<" " 
           <<sz_node_min*(sizeof(SPComplexType)+sizeof(int))/1024.0/1024.0 <<" MB   -   " 
           <<sz_node_max <<" " 
           <<sz_node_max*(sizeof(SPComplexType)+sizeof(int))/1024.0/1024.0 <<" MB \n" 
           <<"    Global: "
           <<tot_sz_global <<" " 
           <<tot_sz_global*(sizeof(SPComplexType)+sizeof(int))/1024.0/1024.0 <<" MB" 
           <<std::endl <<std::endl;

    if(global_reserve)
      reserve_to_fit(Vijkl,sz_global); 
    else
      reserve_to_fit(Vijkl,sz_local); 
  }

  // now calculate fully distributed matrix elements
  for(int nn=0, nb=0, nkcum=0; nn<ngrp; nn++) {

    // just checking
    assert(nkcum==M_split[nn]);
    if(M_split[nn+1]==M_split[nn]) continue;
    int nblk = Qknum[nn];
    long ntermscum=0;
    for( int bi = 0; bi < nblk; bi++, nb++) {
      int nterms = Qksizes[2*nb];      // number of terms in block 
      int nk = Qksizes[2*nb+1];        // number of k-blocks in block
      int k0 = nkcum;                  // first value of k in block
      nkcum+=nk;
      int kN = nkcum;                  // last+1 value
      int NEL0 = (k0<NMO)?NAEA:NAEB;   // number of electrons in this spin block
      assert(nk > 0 && nk <= maxnk );  // just checking

      boost::multi_array_ref<SPComplexType,2> tQk(tQk_shmbuff.data(),extents[nk*NEL0][nvec]);

      Timer_.reset("T0");
      Timer_.start("T0");
      if(sparseQk) {
        if(coreid==0) {
          if(nn == nodeid) {
            auto ka0 = (k0-M_split[nn])*NEL0;
            auto kaN = (k0-M_split[nn]+nk)*NEL0;
            auto n0 = *SpQk.pointers_begin( ka0 );
            auto n1 = *SpQk.pointers_end(kaN);
            int nt_ = static_cast<int>(n1-n0);
            assert(ntermscum==n0);
            assert(nt_==nterms);
            std::copy(std::addressof(*SpQk.non_zero_values_data(ka0)),
                        std::addressof(*SpQk.non_zero_values_data(kaN)),
                        std::addressof(*SptQk.non_zero_values_data()));  
            std::copy(std::addressof(*SpQk.non_zero_indices2_data(ka0)),
                        std::addressof(*SpQk.non_zero_indices2_data(kaN)),
                        std::addressof(*SptQk.non_zero_indices2_data()));  
            for(int i=0, j=ka0; i<=nk*NEL0; i++, j++) {
              SptQk.pointers_begin()[i] = SpQk.pointers_begin()[j]-n0;
              SptQk.pointers_end()[i] = SpQk.pointers_end()[j]-n0;
            }
          }
          TG.Cores().broadcast_value(nterms,nn);  
          TG.Cores().broadcast_n(std::addressof(*SptQk.pointers_begin()),
                                  SptQk.shape()[0],nn);
          TG.Cores().broadcast_n(std::addressof(*SptQk.pointers_end()),
                                  SptQk.shape()[0],nn);
          TG.Cores().broadcast_n(std::addressof(*SptQk.non_zero_values_data()),
                                  nterms,nn);  
          TG.Cores().broadcast_n(std::addressof(*SptQk.non_zero_indices2_data()),
                                  nterms,nn);  
        }
        TG.node_barrier();
        // for safety, keep track of sum
        ntermscum += static_cast<long>(nterms);
      } else {
        if(coreid==0) {
          if(nn == nodeid)
            std::copy(Qk.origin()+bi*maxnk*NEL0*nvec,Qk.origin()+(bi*maxnk+nk)*NEL0*nvec,tQk.origin());
          TG.Cores().broadcast_n(tQk.origin(),nk*NEL0*nvec,nn);
        }
        TG.node_barrier();
      }
      app_log()<<" Loop: " <<nn <<"/" <<ngrp <<" " <<bi <<"/" <<nblk
                 <<" communication: " <<Timer_.total("T0") <<" ";

      boost::multi_array_ref<ComplexType,2> Ta(Ta_shmbuff.data(),extents[nk*NEL0][nrow]);

      Timer_.reset("T0");
      Timer_.start("T0");
      if(type == "SD")
        Qk_x_Rl(walker_type,TG,k0,kN,l0,lN,NMO,NAEA,NAEB,SptQk[{0,std::size_t(nk*NEL0)}],Rl,Ta,Vijkl,cut);
      else if(type == "DD")   
        Qk_x_Rl(walker_type,TG,k0,kN,l0,lN,NMO,NAEA,NAEB,tQk,Rl,Ta,Vijkl,cut);
      Timer_.stop("T0");
      app_log()<<" QxR: " <<Timer_.total("T0") <<std::endl;
//      else if(type == "SS")

    }
  }

  TG.node_barrier();

}

template<class Container = std::vector<std::tuple<int,int,SPComplexType>>>
inline void rotateHijklSymmetric(WALKER_TYPES walker_type, TaskGroup_& TG, Container& Vijkl, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, SpVType_shm_csr_matrix const& V2_fact, const RealType cut, int maximum_buffer_size, bool reserve_to_fit_=true, bool global_reserve = true)
{
  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  bool distribute_Ham  = TG.getNumberOfTGs() > 1;
  if(distribute_Ham)
    APP_ABORT(" Distributed V2_fact not yet implemented. \n");

  int NAEA = Alpha->shape()[0];
  int NMO = Alpha->shape()[1];
  int NAEB = NAEA;
  if( walker_type == COLLINEAR ) NAEB = Beta->shape()[0];

  // <ab||kl> = sum_n Q(a,k,n) * Q(b,l,n) - Q(a,l,n)*Q(b,k,n),
  // where:
  //   Q(a,k,n) = sum_i conj(Amat(i,a)) * V2_fact(ik,n)
  //
  // For type==COLLINEAR, only the alpha/alpha and beta/beta sectors of the
  // double bar integrals are assembled. The alpha/beta coulomb part is assumed   
  // to be calculated from the Choleskt matrix directly in HamiltonianOperations 
  //
  // For parallelization, distribute k over nodes.
  // Build ahead of time Q matrices in shared memory to reduce memory/setup time.
  // Assemble integrals in parallel and fully distributed.
  // Collect on all nodes.
  //    - For distributed hamiltonians, you do not need to worry about keeping contiguous
  //    segments of the hamiltonian. Only that the distribution over nodes is roughly equal.
  //    Write a simple algorithm that balances the number of terms in a TG. 
  //  

  mpi3_SHMBuffer<SPComplexType> tQk_shmbuff(TG.Node(),1);
  mpi3_SHMBuffer<SPComplexType> Qk_shmbuff(TG.Node(),1);

  int NMO2 = (walker_type == CLOSED)?NMO:2*NMO;
  int ngrp = std::min(NMO2,nnodes);
  std::vector<int> M_split(ngrp+1);

  // split the orbitals among processors:  ngrp/2 + ngrp%2 for alpha, ngrp/2 for beta 
  M_split[0]=0;
  if(walker_type == CLOSED) {
    FairDivide(NMO2,ngrp,M_split);
  } else if (walker_type == COLLINEAR) {
    if(ngrp==1)
      APP_ABORT(" Error: Current implementation of rotateHijkl requires at least 2 nodes.\n");
    int nalpha = ngrp/2;
    for(int i=0; i<nalpha; i++) {
      int m0,m1;
      std::tie(m0,m1) = FairDivideBoundary(i,NMO,nalpha);
      M_split[i+1]=m1;
    }
    assert(M_split[nalpha]==NMO);
    for(int i=0, nbeta=ngrp/2+ngrp%2; i<nbeta; i++) {
      int m0,m1;
      std::tie(m0,m1) = FairDivideBoundary(i,NMO,nbeta);
      M_split[i+nalpha+1]=NMO+m1;
    }
    assert(M_split[ngrp]==2*NMO);
  } else if (walker_type == NONCOLLINEAR) {
    APP_ABORT(" Finish. \n");
  }

  // Construct your set of Q(k,a,m), R(l,b,m)
  int l0 = (nodeid<ngrp)?(M_split[nodeid]):(-1);
  int lN = (nodeid<ngrp)?(M_split[nodeid+1]):(-1);
  bool amIAlpha = true; // for simplicity, subset of bands must be either elpha or beta, not mixed
  if( l0 < NMO && (lN-1) < NMO )
    amIAlpha = true;
  else if( l0 >= NMO && lN >= NMO )
    amIAlpha = false;
  else {
    std::cerr<<"l0, lN, nodeid, ngrp, NMO: " <<l0 <<" " <<lN <<" " <<nodeid <<" " <<ngrp <<" " <<NMO <<std::endl;
    std::cerr<<" Error: Current algorithm requires an even number of processors. \n\n\n";
    APP_ABORT(" Error: Current algorithm requires an even number of processors. \n\n\n");
  }
  const int NEL = (amIAlpha)?NAEA:NAEB;

  // create new communicator based on amIAlpha and communicate Qk over it  
  int key = (amIAlpha?1:2);
  if(l0 < 0) key=0;  
  boost::mpi3::communicator comm(TG.Cores().split(key));

  // no cross terms, so alpha/beta teams work independently  
  int norb = lN-l0;
  int maxnorb = 0;
  for(int i=0; i<comm.size(); i++) maxnorb = std::max(maxnorb,M_split[i+1]-M_split[i]);
  int nvec = V2_fact.shape()[1];
  // must gather over heads of TG to get nchol per TG and total # chol vecs 
  // Rl(k, a, m), k:[0:NMO2], a:[0:NAEA], m:[0:nvec] 

  app_log()<<" Approximate memory usage for half-rotated Hamiltonian construction: \n"
           <<"   max. number of orbital in a node: " <<maxnorb <<"\n"
           <<"   Q matrix size (assuming dense) = maxnorb * nup * ncholvecs complex numbers = "
           <<sizeof(SPComplexType)*maxnorb*NEL*nvec/1024.0/1024.0 <<" MB \n"
           <<"   Maximum size of communication buffer: " <<maximum_buffer_size <<" MB" <<std::endl;

  const int nrow = norb * NEL; 
  const int ncol = nvec;
  int dummy_nrow=nrow, dummy_ncol=ncol;
  int mat_size = nrow*ncol; 
  if(nodeid < ngrp) {
    Qk_shmbuff.resize(mat_size);
    if(coreid==0) std::fill_n(Qk_shmbuff.data(),Qk_shmbuff.size(),SPComplexType(0.0));
  } else {
    dummy_nrow=dummy_ncol=0;
  }
  TG.node_barrier();

  boost::multi_array_ref<SPComplexType,2> Qk(Qk_shmbuff.data(),extents[dummy_nrow][dummy_ncol]);  
  dummy_nrow=nrow; dummy_ncol=ncol;

  if(distribute_Ham) {
   APP_ABORT(" Finish THIS (43)!!! \n\n\n");
  } else {

    //   Q(k,a,n) = sum_i conj(Amat(i,a)) * V2_fact(ik,n)
    ma_rotate::halfRotateCholeskyMatrix(walker_type,TG,l0,lN,Qk,Alpha,Beta,V2_fact,false,false,cut);

  }

  TG.node_barrier();
  // let the maximum message be maximum_buffer_size MB
  int maxnt = std::max(1,static_cast<int>(std::floor(maximum_buffer_size*1024.0*1024.0/sizeof(SPComplexType))));

  // control the size of the MPI messages.
  // communicate blocks of Qk up to maxnt terms  
  std::vector<int> nkbounds;          // local bounds for communication  
  std::vector<int> Qknum(comm.size());     // number of blocks per node
  std::vector<int> Qksizes;           // number of terms and number of k vectors in block for all nodes
  if(coreid==0) {
    int n_=0, ntcnt=0, n0=0;
    for(int i=0; i<norb; i++) {
      int ntt;
      ntt = nvec*NEL;
      assert(ntt < maxnt);
      if(ntcnt+ntt > maxnt) {
        nkbounds.push_back(ntcnt);
        nkbounds.push_back(i-n0);
        ntcnt=ntt;
        n0=i;
        n_++;
      } else {
        ntcnt+=ntt;
      }
    }
    if(ntcnt > 0) {
      // push last block
      n_++;
      nkbounds.push_back(ntcnt);
      nkbounds.push_back(norb-n0);
    }

    MPI_Allgather(&n_,1,MPI_INT,Qknum.data(),1,MPI_INT,comm.impl_);

    int ntt = std::accumulate(Qknum.begin(),Qknum.end(),0);
    Qksizes.resize(2*ntt);

    std::vector<int> cnts(comm.size());
    std::vector<int> disp(comm.size());
    int cnt=0;
    for(int i=0; i<comm.size(); i++) {
      cnts[i] = Qknum[i]*2;
      disp[i]=cnt;
      cnt+=cnts[i];
    }
    MPI_Allgatherv(nkbounds.data(),nkbounds.size(),MPI_INT,Qksizes.data(),cnts.data(),disp.data(),MPI_INT,comm.impl_ );

  }

  MPI_Bcast(Qknum.data(),comm.size(),MPI_INT,0,TG.Node().impl_);
  int ntt = std::accumulate(Qknum.begin(),Qknum.end(),0);
  if(coreid!=0)
    Qksizes.resize(2*ntt);
  MPI_Bcast(Qksizes.data(),Qksizes.size(),MPI_INT,0,TG.Node().impl_);

// store {nterms,nk} for all nodes 
// use it to know communication pattern 

  int maxnk = 0;        // maximum number of k vectors in communication block 
  long maxqksize = 0;   // maximum size of communication block
  for(int i=0; i<ntt; i++) {
    if(Qksizes[2*i] > maxqksize) maxqksize = Qksizes[2*i];
    if(Qksizes[2*i+1] > maxnk) maxnk = Qksizes[2*i+1];
  }

  app_log()<<"   Maximum number of k blocks communicated in an iteration: " <<maxnk <<"\n"
           <<"   Temporary integral matrix Ta: " <<norb*NEL*maxnk*NEL*sizeof(SPComplexType)/1024.0/1024.0 <<" MB " <<std::endl;

  // temporary shared memory space for local "dense" result
  mpi3_SHMBuffer<SPComplexType> Ta_shmbuff(TG.Node(),norb*NEL*maxnk*NEL);
  tQk_shmbuff.resize(maxnk * NEL * nvec);

  myTimer Timer_;

  if(reserve_to_fit_) {
    // count and resize container
    std::vector<std::size_t> sz_local;
    if(walker_type==CLOSED) sz_local.resize(NMO*NAEA);
    else if(walker_type==COLLINEAR) sz_local.resize(NMO*(NAEA+NAEB));
    else if(walker_type==NONCOLLINEAR) sz_local.resize(2*NMO*(NAEA+NAEB));
    
    int nkcum=(amIAlpha?0:NMO);
    for(int nn=0, nb=0; nn<comm.size(); nn++) {

      // just checking
      assert(nkcum==M_split[nn]);
      if(M_split[nn+1]==M_split[nn]) continue;
      int nblk = Qknum[nn];
      long ntermscum=0;
      for( int bi = 0; bi < nblk; bi++, nb++) {
        int nterms = Qksizes[2*nb];      // number of terms in block 
        int nk = Qksizes[2*nb+1];        // number of k-blocks in block
        int k0 = nkcum;                  // first value of k in block
        nkcum+=nk;
        int kN = nkcum;                  // last+1 value
        assert(nk > 0 && nk <= maxnk );  // just checking

        Timer_.reset("T0");
        Timer_.start("T0");
        boost::multi_array_ref<SPComplexType,2> tQk(tQk_shmbuff.data(),extents[nk*NEL][nvec]);
        if(coreid==0) {
          if(nn == comm.rank())
            std::copy(Qk.origin()+bi*maxnk*NEL*nvec,Qk.origin()+(bi*maxnk+nk)*NEL*nvec,tQk.origin());
          comm.broadcast_n(tQk.origin(),nk*NEL*nvec,nn);
        }
        TG.node_barrier();
        app_log()<<" Loop: " <<nn <<"/" <<comm.size() <<" " <<bi <<"/" <<nblk
                 <<" communication: " <<Timer_.total("T0") <<" ";

        boost::multi_array_ref<ComplexType,2> Ta(Ta_shmbuff.data(),extents[nk*NEL][nrow]);

        Timer_.reset("T0");
        Timer_.start("T0");
        count_Qk_x_Qk_symmetric(walker_type,TG,sz_local,k0,kN,l0,lN,NMO,NAEA,NAEB,tQk,Qk,Ta,cut);
        Timer_.stop("T0");
        app_log()<<" QxR: " <<Timer_.total("T0") <<std::endl;

      }
    }

    std::size_t tot_sz_local = std::accumulate(sz_local.begin(),sz_local.end(),std::size_t(0));

    std::vector<std::size_t> sz_global(sz_local.size());
    TG.Global().all_reduce_n(sz_local.begin(),sz_local.size(),sz_global.begin(),std::plus<>());
    std::size_t tot_sz_global = std::accumulate(sz_global.begin(),sz_global.end(),std::size_t(0));

    std::vector<std::size_t> sz_node(sz_local.size());
    TG.Node().all_reduce_n(sz_local.begin(),sz_local.size(),sz_node.begin(),std::plus<>());
    std::size_t tot_sz_node = std::accumulate(sz_node.begin(),sz_node.end(),std::size_t(0));

    std::size_t sz_node_min= TG.Global().all_reduce_value(tot_sz_node,boost::mpi3::min<>());
    std::size_t sz_node_max= TG.Global().all_reduce_value(tot_sz_node,boost::mpi3::max<>());
    std::size_t sz_local_min= TG.Global().all_reduce_value(tot_sz_local,boost::mpi3::min<>());
    std::size_t sz_local_max= TG.Global().all_reduce_value(tot_sz_local,boost::mpi3::max<>());

    app_log()<<"  Number of terms in Vijkl: \n" 
           <<"    Local: (min/max)" 
           <<sz_local_min <<" " 
           <<sz_local_min*(sizeof(SPComplexType)+sizeof(int))/1024.0/1024.0 <<" MB  -  " 
           <<sz_local_max <<" " 
           <<sz_local_max*(sizeof(SPComplexType)+sizeof(int))/1024.0/1024.0 <<" MB \n" 
           <<"    Node (min/max): " 
           <<sz_node_min <<" " 
           <<sz_node_min*(sizeof(SPComplexType)+sizeof(int))/1024.0/1024.0 <<" MB   -   " 
           <<sz_node_max <<" " 
           <<sz_node_max*(sizeof(SPComplexType)+sizeof(int))/1024.0/1024.0 <<" MB \n" 
           <<"    Global: "
           <<tot_sz_global <<" " 
           <<tot_sz_global*(sizeof(SPComplexType)+sizeof(int))/1024.0/1024.0 <<" MB" 
           <<std::endl <<std::endl;

    if(global_reserve)
      reserve_to_fit(Vijkl,sz_global); 
    else
      reserve_to_fit(Vijkl,sz_local); 
  }

  // now calculate fully distributed matrix elements
  int nkcum=(amIAlpha?0:NMO);
  for(int nn=0, nb=0; nn<comm.size(); nn++) {

    // just checking
    assert(nkcum==M_split[nn]);
    if(M_split[nn+1]==M_split[nn]) continue;
    int nblk = Qknum[nn];
    long ntermscum=0;
    for( int bi = 0; bi < nblk; bi++, nb++) {
      int nterms = Qksizes[2*nb];      // number of terms in block 
      int nk = Qksizes[2*nb+1];        // number of k-blocks in block
      int k0 = nkcum;                  // first value of k in block
      nkcum+=nk;
      int kN = nkcum;                  // last+1 value
      assert(nk > 0 && nk <= maxnk );  // just checking

      boost::multi_array_ref<SPComplexType,2> tQk(tQk_shmbuff.data(),extents[nk*NEL][nvec]);

      Timer_.reset("T0");
      Timer_.start("T0");
      if(coreid==0) {
        if(nn == nodeid)
          std::copy(Qk.origin()+bi*maxnk*NEL*nvec,Qk.origin()+(bi*maxnk+nk)*NEL*nvec,tQk.origin());
        comm.broadcast_n(tQk.origin(),nk*NEL*nvec,nn);
      }
      TG.node_barrier();
      app_log()<<" Loop: " <<nn <<"/" <<comm.size() <<" " <<bi <<"/" <<nblk
                 <<" communication: " <<Timer_.total("T0") <<" ";

      boost::multi_array_ref<ComplexType,2> Ta(Ta_shmbuff.data(),extents[nk*NEL][nrow]);

      Timer_.reset("T0");
      Timer_.start("T0");
      Qk_x_Qk_symmetric(walker_type,TG,k0,kN,l0,lN,NMO,NAEA,NAEB,tQk,Qk,Ta,Vijkl,cut);
      Timer_.stop("T0");
      app_log()<<" QxR: " <<Timer_.total("T0") <<std::endl;

    }
  }

  TG.node_barrier();

}

}

}

#endif

