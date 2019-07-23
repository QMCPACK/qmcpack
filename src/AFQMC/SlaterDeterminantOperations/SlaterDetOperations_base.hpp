//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_BASE_HPP
#define QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_BASE_HPP

#include<fstream>

#include "AFQMC/config.h"
#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/csr_blas.hpp"
#include "AFQMC/SlaterDeterminantOperations/mixed_density_matrix.hpp"
#include "AFQMC/SlaterDeterminantOperations/apply_expM.hpp"

#include "type_traits/scalar_traits.h"

namespace qmcplusplus
{

namespace afqmc
{

template<class AllocType> 
class SlaterDetOperations_base 
{
  public:

    using Alloc = AllocType;
    using T = typename Alloc::value_type;
    using R = typename remove_complex<T>::value_type; 
    using pointer = typename Alloc::pointer;
    using const_pointer = typename Alloc::const_pointer;
    using IAlloc = typename Alloc::template rebind<int>::other;
    using RAlloc = typename Alloc::template rebind<R>::other;
    using Rpointer = typename RAlloc::pointer;

    using IVector = boost::multi::array<int,1,IAlloc>;  
    using TVector = boost::multi::array<T,1,Alloc>;  
    using RVector = boost::multi::array<R,1,RAlloc>;  
    using TVector_ref = boost::multi::array_ref<T,1,pointer>;  
    using TMatrix = boost::multi::array<T,2,Alloc>;  
    using TMatrix_ref = boost::multi::array_ref<T,2,pointer>;  
    using TTensor_ref = boost::multi::array_ref<T,3,pointer>;  

    SlaterDetOperations_base(Alloc alloc_={}):
      allocator_(alloc_),
      iallocator_(alloc_),
      WORK(iextensions<1u>{0},allocator_),
      IWORK(iextensions<1u>{0},iallocator_),
      TAU(iextensions<1u>{0},allocator_),
      TBuff(iextensions<1u>{0},allocator_)
    {
    }

    SlaterDetOperations_base(int NMO, int NAEA, Alloc alloc_={}):
      allocator_(alloc_),
      iallocator_(alloc_),
      rallocator_(alloc_),
      WORK(iextensions<1u>{0},allocator_),  
      IWORK(iextensions<1u>{NMO+1},iallocator_),  
      RWORK(iextensions<1u>{NMO+1},rallocator_),  
      TAU(iextensions<1u>{NMO},allocator_),  
      TBuff(iextensions<1u>{NMO*NMO},allocator_)
    {

      // only used to determine size of WORK  
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNM(TBuff.data(), {NAEA,NMO});
      TMatrix_ref TMN(TBuff.data(), {NMO,NAEA});

      // reserve enough space in lapack's work array
      // Make sure it is large enough for:
      // 1. getri( TNN )
      //  2. geqrf( TNM )
      int mem_needs = std::max(ma::getri_optimal_workspace_size(TNN),
                                  ma::geqrf_optimal_workspace_size(TNM));  
      //  3. gqr( TNM )
      mem_needs = std::max(mem_needs, ma::gqr_optimal_workspace_size(TNM) );
      //  4. gelqf( TMN )
      mem_needs = std::max(mem_needs, ma::gelqf_optimal_workspace_size(TMN) ); 
      //  5. glq( TMN )
      mem_needs = std::max(mem_needs, ma::glq_optimal_workspace_size(TMN) ); 
      //  6. trf( TNN )
      mem_needs = std::max(mem_needs, ma::getrf_optimal_workspace_size(TNN) ); 
      //  7. svd( TNN )
      mem_needs = std::max(mem_needs, ma::gesvd_optimal_workspace_size(TNN) ); 
      mem_needs = std::max(mem_needs, NMO); 
      WORK.reextent( iextensions<1u>{mem_needs} );   
    }

    ~SlaterDetOperations_base() {}

    SlaterDetOperations_base(const SlaterDetOperations_base& other) = delete;
    SlaterDetOperations_base(SlaterDetOperations_base&& other) = default;
    SlaterDetOperations_base& operator=(const SlaterDetOperations_base& other) = delete;
    SlaterDetOperations_base& operator=(SlaterDetOperations_base&& other) = default;

    template<class MatA, class MatB, class MatC>
    T MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, T LogOverlapFactor, bool compact, bool herm=true) {
      int NMO = (herm?hermA.size(1):hermA.size(0));
      int NAEA = (herm?hermA.size(0):hermA.size(1));
      set_buffer(NAEA*NAEA + NAEA*NMO);
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNM(TBuff.data()+TNN.num_elements(), {NAEA,NMO});
      return SlaterDeterminantOperations::base::MixedDensityMatrix<T>(hermA,B,std::forward<MatC>(C),LogOverlapFactor,TNN,TNM,IWORK,WORK,compact,herm);
    }

    template<class MatA, class MatC>
    T MixedDensityMatrix(const MatA& A, MatC&& C, T LogOverlapFactor, bool compact=false) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      set_buffer(NAEA*NAEA + NAEA*NMO);
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNM(TBuff.data()+TNN.num_elements(), {NAEA,NMO});
      return SlaterDeterminantOperations::base::MixedDensityMatrix<T>(A,A,std::forward<MatC>(C),LogOverlapFactor,TNN,TNM,IWORK,WORK,compact,false);
    }

    template<class MatA, class MatB, class MatC>
    T MixedDensityMatrix_noHerm(const MatA& A, const MatB& B, MatC&& C, T LogOverlapFactor, bool compact=false, bool useSVD = false) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      if(useSVD) {
        if( RWORK.num_elements() < 6*NAEA+1 ) RWORK.reextent(iextensions<1u>{6*NAEA+1});
        set_buffer(2*NAEA*NAEA + 2*NAEA*NMO);
        TMatrix_ref TNN1(TBuff.data(), {NAEA,NAEA});
        TMatrix_ref TNN2(TNN1.origin()+TNN1.num_elements(), {NAEA,NAEA});
        TMatrix_ref TNM(TNN2.origin()+TNN2.num_elements(), {NAEA,NMO});
        TMatrix_ref TMN(TNM.origin()+TNM.num_elements(), {NMO,NAEA});
        return SlaterDeterminantOperations::base::MixedDensityMatrix_noHerm_wSVD<T>(A,B,std::forward<MatC>(C),LogOverlapFactor,RWORK,TNN1,TNN2,TMN,TNM,IWORK,WORK,compact);
      } else {
        set_buffer(NAEA*NAEA + NAEA*NMO);
        TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
        TMatrix_ref TNM(TBuff.data()+TNN.num_elements(), {NAEA,NMO});
        return SlaterDeterminantOperations::base::MixedDensityMatrix_noHerm<T>(A,B,std::forward<MatC>(C),LogOverlapFactor,TNN,TNM,IWORK,WORK,compact);
      }
    }

    template<class integer, class MatA, class MatB, class MatC, class MatQ>
    T MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C, T LogOverlapFactor,  
                                    integer* ref, MatQ&& QQ0, bool compact=false) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      int NMO = B.size(0);
      assert(hermA.size(1)==B.size(0));
      assert(QQ0.size(0)==Nact);
      assert(QQ0.size(1)==NEL);
      set_buffer(NEL*NEL + Nact*NEL + NMO*NEL);
      TMatrix_ref TNN(TBuff.data(), {NEL,NEL});
      TMatrix_ref TAB(TNN.data()+TNN.num_elements(), {Nact,NEL});
      TMatrix_ref TNM(TAB.data()+TAB.num_elements(), {NEL,NMO});
      return SlaterDeterminantOperations::base::MixedDensityMatrixForWoodbury<T>(hermA,B,std::forward<MatC>(C),LogOverlapFactor,std::forward<MatQ>(QQ0),ref,TNN,TAB,TNM,IWORK,WORK,compact);
    }

    template<class integer, class MatA, class MatB, class MatC>
    T MixedDensityMatrixFromConfiguration(const MatA& hermA, const MatB& B, MatC&& C, T LogOverlapFactor, 
                                    integer* ref, bool compact=false) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      int NMO = B.size(0);
      assert(hermA.size(1)==B.size(0));
      set_buffer(NEL*NEL + Nact*NEL + NMO*NEL);
      TMatrix_ref TNN(TBuff.data(), {NEL,NEL});
      TMatrix_ref TAB(TNN.data()+TNN.num_elements(), {Nact,NEL});
      TMatrix_ref TNM(TAB.data()+TAB.num_elements(), {NEL,NMO});
      return SlaterDeterminantOperations::base::MixedDensityMatrixFromConfiguration<T>(hermA,B,std::forward<MatC>(C),LogOverlapFactor,ref,TNN,TAB,TNM,IWORK,WORK,compact);
    }

    template<class MatA>
    T Overlap(const MatA& A, T LogOverlapFactor) {
      int NAEA = A.size(1);
      set_buffer(2*NAEA*NAEA);
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNN2(TBuff.data()+TNN.num_elements(), {NAEA,NAEA});
      return SlaterDeterminantOperations::base::Overlap<T>(A,A,LogOverlapFactor,TNN,IWORK,TNN2,false);
    }

    template<class MatA, class MatB>
    T Overlap(const MatA& hermA, const MatB& B, T LogOverlapFactor, bool herm=true) {
      int NAEA = (herm?hermA.size(0):hermA.size(1));
      set_buffer(2*NAEA*NAEA);
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNN2(TBuff.data()+TNN.num_elements(), {NAEA,NAEA});
      return SlaterDeterminantOperations::base::Overlap<T>(hermA,B,LogOverlapFactor,TNN,IWORK,TNN2,herm);
    } 

    template<class MatA, class MatB>
    T Overlap_noHerm(const MatA& A, const MatB& B, T LogOverlapFactor) {
      int NAEA = A.size(1);
      set_buffer(2*NAEA*NAEA);
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNN2(TBuff.data()+TNN.num_elements(), {NAEA,NAEA});
      return SlaterDeterminantOperations::base::Overlap<T>(A,B,LogOverlapFactor,TNN,IWORK,TNN2,false);
    }

    // routines for PHMSD
    template<typename integer, class MatA, class MatB, class MatC>
    T OverlapForWoodbury(const MatA& hermA, const MatB& B, T LogOverlapFactor, integer* ref, MatC&& QQ0) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      assert(hermA.size(1)==B.size(0));
      assert(QQ0.size(0)==Nact);  
      assert(QQ0.size(1)==NEL);  
      set_buffer(NEL*NEL + Nact+NEL);
      TMatrix_ref TNN(TBuff.data(), {NEL,NEL});
      TMatrix_ref TMN(TBuff.data()+TNN.num_elements(), {Nact,NEL});
      return SlaterDeterminantOperations::base::OverlapForWoodbury<T>(hermA,B,LogOverlapFactor,std::forward<MatC>(QQ0),ref,TNN,TMN,IWORK,WORK);
    }    

    template<class Mat, class MatP1, class MatV>
    void Propagate(Mat&& A, const MatP1& P1, const MatV& V, int order=6, char TA='N') {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      set_buffer(3*NMO*NAEA);
      TMatrix_ref TMN(TBuff.data(), {NMO,NAEA});
      TMatrix_ref T1(TMN.data()+TMN.num_elements(), {NMO,NAEA});
      TMatrix_ref T2(T1.data()+T1.num_elements(), {NMO,NAEA});
      using ma::T;
      using ma::H;
      if(TA=='H' || TA=='h') {
        ma::product(ma::H(P1),std::forward<Mat>(A),TMN);
        SlaterDeterminantOperations::base::apply_expM(V,TMN,T1,T2,order,TA);
        ma::product(ma::H(P1),TMN,std::forward<Mat>(A));
      } else if(TA=='T' || TA=='t') {
        ma::product(ma::T(P1),std::forward<Mat>(A),TMN);
        SlaterDeterminantOperations::base::apply_expM(V,TMN,T1,T2,order,TA);
        ma::product(ma::T(P1),TMN,std::forward<Mat>(A));
      } else {
        ma::product(P1,std::forward<Mat>(A),TMN);
        SlaterDeterminantOperations::base::apply_expM(V,TMN,T1,T2,order);
        ma::product(P1,TMN,std::forward<Mat>(A));
      }  
    }

// NOTE: Move to a single buffer for all TNXs, to be able to reuse buffers 
//       more efficiently and to eventually share them with HamOps  

    // need to check if this is equivalent to QR!!!
    template<class Mat>
    T Orthogonalize(Mat&& A, T LogOverlapFactor) {
#ifdef ENABLE_CUDA
      // QR on the transpose
      int NMO = A.size(0);
      int NAEA = A.size(1);
      set_buffer(NMO*NAEA + NMO);
      TMatrix_ref AT(TBuff.origin(), {NAEA,NMO});
      TVector_ref scl(AT.origin() + AT.num_elements(), {NMO});
      ma::transpose(A,AT);   
      ma::geqrf(AT,TAU,WORK);
      using ma::determinant_from_geqrf;
      using ma::scale_columns;  
      T res = determinant_from_geqrf(AT.size(0),AT.origin(),AT.stride(0),scl.origin(),LogOverlapFactor);
      ma::gqr(AT,TAU,WORK);
      ma::transpose(AT,A);   
      scale_columns(A.size(0),A.size(1),A.origin(),A.stride(0),scl.origin());
#else
      ma::gelqf(std::forward<Mat>(A),TAU,WORK);
      T res(0.0); 
       for (int i = 0; i < A.size(1); i++) { 
        if (real(A[i][i]) < 0) 
          IWORK[i]=-1; 
        else 
          IWORK[i]=1; 
        res += std::log(T(IWORK[i])*A[i][i]);
      }
      res = std::exp(res - LogOverlapFactor);
      ma::glq(std::forward<Mat>(A),TAU,WORK);
      for(int i=0; i<A.size(0); ++i)
        for(int j=0; j<A.size(1); ++j)
          A[i][j] *= T(IWORK[j]);
#endif
      return res;
    }

  protected:

    Alloc allocator_;
    IAlloc iallocator_;
    RAlloc rallocator_;

    TVector WORK;
    IVector IWORK;
    RVector RWORK;

    // Vector used in QR routines 
    TVector TAU;

    TVector TBuff;

    void set_buffer(size_t N) {
      if(TBuff.num_elements() < N) 
        TBuff = std::move(TVector(iextensions<1u>{N}));
      using std::fill_n;
      fill_n(TBuff.origin(),N,T(0.0));
    }

};

}

}

#endif
