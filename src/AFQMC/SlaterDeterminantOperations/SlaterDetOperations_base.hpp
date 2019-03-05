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
    using pointer = typename Alloc::pointer;
    using const_pointer = typename Alloc::const_pointer;
    using IAlloc = typename Alloc::template rebind<int>::other;;

    using IVector = boost::multi::array<int,1,IAlloc>;  
    using TVector = boost::multi::array<T,1,Alloc>;  
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
      WORK(iextensions<1u>{0},allocator_),  
      IWORK(iextensions<1u>{NMO+1},iallocator_),  
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
      WORK.reextent( iextensions<1u>{mem_needs} );   
    }

    ~SlaterDetOperations_base() {}

    SlaterDetOperations_base(const SlaterDetOperations_base& other) = delete;
    SlaterDetOperations_base(SlaterDetOperations_base&& other) = default;
    SlaterDetOperations_base& operator=(const SlaterDetOperations_base& other) = delete;
    SlaterDetOperations_base& operator=(SlaterDetOperations_base&& other) = default;

    template<class MatA, class MatB, class MatC>
    void MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, T* res, bool compact=false) {
      int NMO = hermA.size(1);
      int NAEA = hermA.size(0);
      set_buffer(NAEA*NAEA + NAEA*NMO);
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNM(TBuff.data()+TNN.num_elements(), {NAEA,NMO});
      SlaterDeterminantOperations::base::MixedDensityMatrix<T>(hermA,B,std::forward<MatC>(C),res,TNN,TNM,IWORK,WORK,compact);
    }

    template<class MatA, class MatC>
    void MixedDensityMatrix(const MatA& A, MatC&& C, T* res, bool compact=false) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      set_buffer(NAEA*NAEA + NAEA*NMO);
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNM(TBuff.data()+TNN.num_elements(), {NAEA,NMO});
      SlaterDeterminantOperations::base::MixedDensityMatrix_noHerm<T>(A,A,std::forward<MatC>(C),res,TNN,TNM,IWORK,WORK,compact);
    }

    template<class MatA, class MatB, class MatC>
    void MixedDensityMatrix_noHerm(const MatA& A, const MatB& B, MatC&& C, T* res, bool compact=false) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      set_buffer(NAEA*NAEA + NAEA*NMO);
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNM(TBuff.data()+TNN.num_elements(), {NAEA,NMO});
      SlaterDeterminantOperations::base::MixedDensityMatrix_noHerm<T>(A,B,std::forward<MatC>(C),res,TNN,TNM,IWORK,WORK,compact);
    }

    template<class integer, class MatA, class MatB, class MatC, class MatQ>
    void MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C, T* res, 
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
      SlaterDeterminantOperations::base::MixedDensityMatrixForWoodbury<T>(hermA,B,std::forward<MatC>(C),res,std::forward<MatQ>(QQ0),ref,TNN,TAB,TNM,IWORK,WORK,compact);
    }

    template<class integer, class MatA, class MatB, class MatC>
    void MixedDensityMatrixFromConfiguration(const MatA& hermA, const MatB& B, MatC&& C, T* res,
                                    integer* ref, bool compact=false) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      int NMO = B.size(0);
      assert(hermA.size(1)==B.size(0));
      set_buffer(NEL*NEL + Nact*NEL + NMO*NEL);
      TMatrix_ref TNN(TBuff.data(), {NEL,NEL});
      TMatrix_ref TAB(TNN.data()+TNN.num_elements(), {Nact,NEL});
      TMatrix_ref TNM(TAB.data()+TAB.num_elements(), {NEL,NMO});
      SlaterDeterminantOperations::base::MixedDensityMatrixFromConfiguration<T>(hermA,B,std::forward<MatC>(C),res,ref,TNN,TAB,TNM,IWORK,WORK,compact);
    }

    template<class MatA>
    void Overlap(const MatA& A, T* res) {
      int NAEA = A.size(1);
      set_buffer(2*NAEA*NAEA);
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNN2(TBuff.data()+TNN.num_elements(), {NAEA,NAEA});
      SlaterDeterminantOperations::base::Overlap_noHerm<T>(A,A,res,TNN,IWORK,TNN2);
    }

    template<class MatA, class MatB>
    void Overlap(const MatA& hermA, const MatB& B, T* res) {
      int NAEA = hermA.size(0);
      set_buffer(2*NAEA*NAEA);
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNN2(TBuff.data()+TNN.num_elements(), {NAEA,NAEA});
      SlaterDeterminantOperations::base::Overlap<T>(hermA,B,res,TNN,IWORK,TNN2);
    } 

    template<class MatA, class MatB>
    void Overlap_noHerm(const MatA& A, const MatB& B, T* res) {
      int NAEA = A.size(1);
      set_buffer(2*NAEA*NAEA);
      TMatrix_ref TNN(TBuff.data(), {NAEA,NAEA});
      TMatrix_ref TNN2(TBuff.data()+TNN.num_elements(), {NAEA,NAEA});
      SlaterDeterminantOperations::base::Overlap_noHerm<T>(A,B,res,TNN,IWORK,TNN2);
    }

    // routines for PHMSD
    template<typename integer, class MatA, class MatB, class MatC>
    void OverlapForWoodbury(const MatA& hermA, const MatB& B, T* res, integer* ref, MatC&& QQ0) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      assert(hermA.size(1)==B.size(0));
      assert(QQ0.size(0)==Nact);  
      assert(QQ0.size(1)==NEL);  
      set_buffer(NEL*NEL + Nact+NEL);
      TMatrix_ref TNN(TBuff.data(), {NEL,NEL});
      TMatrix_ref TMN(TBuff.data()+TNN.num_elements(), {Nact,NEL});
      SlaterDeterminantOperations::base::OverlapForWoodbury<T>(hermA,B,res,std::forward<MatC>(QQ0),ref,TNN,TMN,IWORK,WORK);
    }    

    template<class Mat, class MatP1, class MatV>
    void Propagate(Mat&& A, const MatP1& P1, const MatV& V, int order=6) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      set_buffer(3*NMO*NAEA);
      TMatrix_ref TMN(TBuff.data(), {NMO,NAEA});
      TMatrix_ref T1(TMN.data()+TMN.num_elements(), {NMO,NAEA});
      TMatrix_ref T2(T1.data()+T1.num_elements(), {NMO,NAEA});
      ma::product(P1,std::forward<Mat>(A),TMN);
      SlaterDeterminantOperations::base::apply_expM(V,TMN,T1,T2,order);
      ma::product(P1,TMN,std::forward<Mat>(A));
    }

// NOTE: Move to a single buffer for all TNXs, to be able to reuse buffers 
//       more efficiently and to eventually share them with HamOps  

    // need to check if this is equivalent to QR!!!
    template<class Mat>
    void Orthogonalize(Mat&& A, T* res=nullptr) {
#ifdef QMC_CUDA
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
      if(res != nullptr)
        determinant_from_geqrf(AT.size(0),AT.origin(),AT.stride(0),scl.origin(),res);
      else {
        determinant_from_geqrf(AT.size(0),AT.origin(),AT.stride(0),scl.origin());
      }
      ma::gqr(AT,TAU,WORK);
      ma::transpose(AT,A);   
      scale_columns(A.size(0),A.size(1),A.origin(),A.stride(0),scl.origin());
#else
      ma::gelqf(std::forward<Mat>(A),TAU,WORK);
      if(res != nullptr) {
        *res = T(1.0); 
        for (int i = 0; i < A.size(1); i++) { 
          if (real(A[i][i]) < 0) 
            IWORK[i]=-1; 
          else 
            IWORK[i]=1; 
          *res *= T(IWORK[i])*A[i][i];
        }
      } else {
        for (int i = 0; i < A.size(1); i++) { 
          if (real(A[i][i]) < 0) 
            IWORK[i]=-1; 
          else 
            IWORK[i]=1; 
        }
      }
      ma::glq(std::forward<Mat>(A),TAU,WORK);
      for(int i=0; i<A.size(0); ++i)
        for(int j=0; j<A.size(1); ++j)
          A[i][j] *= T(IWORK[j]);
#endif
    }

  protected:

    Alloc allocator_;
    IAlloc iallocator_;

    TVector WORK;
    IVector IWORK;

    // Vector used in QR routines 
    TVector TAU;

    TVector TBuff;

    void set_buffer(size_t N) {
      if(TBuff.num_elements() < N)
        TBuff = std::move(TVector(iextensions<1u>{N}));
    }

};

}

}

#endif
