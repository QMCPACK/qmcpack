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
    using TMatrix = boost::multi::array<T,2,Alloc>;  
    using TMatrix_ref = boost::multi::array_ref<T,2,pointer>;  

    SlaterDetOperations_base(Alloc alloc_={}):
      allocator_(alloc_),
      iallocator_(alloc_),
      WORK(extensions<1u>{0},allocator_),
      IWORK(extensions<1u>{0},iallocator_),
      TAU(extensions<1u>{0},allocator_),
      TMat_NN(extensions<2u>{0,0},allocator_),
      TMat_NM(extensions<2u>{0,0},allocator_),
      TMat_MN(extensions<2u>{0,0},allocator_),
      TMat_MM(extensions<2u>{0,0},allocator_),
      TMat_MM2(extensions<2u>{0,0},allocator_),
      TMat_MM3(extensions<2u>{0,0},allocator_)
    {
    }

    SlaterDetOperations_base(int NMO, int NAEA, Alloc alloc_={}):
      allocator_(alloc_),
      iallocator_(alloc_),
      WORK(extensions<1u>{0},allocator_),  
      IWORK(extensions<1u>{NMO+1},iallocator_),  
      TAU(extensions<1u>{NMO},allocator_),  
      TMat_NN(extensions<2u>{NAEA,NAEA},allocator_),
      TMat_NM(extensions<2u>{NAEA,NMO},allocator_),
      TMat_MN(extensions<2u>{NMO,NAEA},allocator_),
      TMat_MM(extensions<2u>{NMO,NMO},allocator_),
      TMat_MM2(extensions<2u>{NMO,NMO},allocator_),
      TMat_MM3(extensions<2u>{NMO,NMO},allocator_)
    {
      // reserve enough space in lapack's work array
      // Make sure it is large enough for:
      // 1. getri( TMat_NN )
      //  2. geqrf( TMat_NM )
      int mem_needs = std::max(ma::getri_optimal_workspace_size(TMat_NN),
                                  ma::geqrf_optimal_workspace_size(TMat_NM));  
      //  3. gqr( TMat_NM )
      mem_needs = std::max(mem_needs, ma::gqr_optimal_workspace_size(TMat_NM) );
      //  4. gelqf( TMat_MN )
      mem_needs = std::max(mem_needs, ma::gelqf_optimal_workspace_size(TMat_MN) ); 
      //  5. glq( TMat_MN )
      mem_needs = std::max(mem_needs, ma::glq_optimal_workspace_size(TMat_MN) ); 
      WORK.reextent( extensions<1u>{mem_needs} );   
    }

    ~SlaterDetOperations_base() {}

    SlaterDetOperations_base(const SlaterDetOperations_base& other) = delete;
    SlaterDetOperations_base(SlaterDetOperations_base&& other) = default;
    SlaterDetOperations_base& operator=(const SlaterDetOperations_base& other) = delete;
    SlaterDetOperations_base& operator=(SlaterDetOperations_base&& other) = default;

    template<class MatA, class MatB, class MatC>
    void MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, T* res, bool compact=false) {
      int NMO = hermA.shape()[1];
      int NAEA = hermA.shape()[0];
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      TMatrix_ref TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NMO);
      TMatrix_ref TNM(TMat_NM.data(), {NAEA,NMO});
      SlaterDeterminantOperations::base::MixedDensityMatrix<T>(hermA,B,std::forward<MatC>(C),res,TNN,TNM,IWORK,WORK,compact);
    }

    template<class MatA, class MatC>
    void MixedDensityMatrix(const MatA& A, MatC&& C, T* res, bool compact=false) {
      int NMO = A.shape()[0];
      int NAEA = A.shape()[1];
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      TMatrix_ref TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NMO);
      TMatrix_ref TNM(TMat_NM.data(), {NAEA,NMO});
      SlaterDeterminantOperations::base::MixedDensityMatrix_noHerm<T>(A,A,std::forward<MatC>(C),res,TNN,TNM,IWORK,WORK,compact);
    }

    template<class MatA, class MatB, class MatC>
    void MixedDensityMatrix_noHerm(const MatA& A, const MatB& B, MatC&& C, T* res, bool compact=false) {
      int NMO = A.shape()[0];
      int NAEA = A.shape()[1];
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      TMatrix_ref TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NMO);
      TMatrix_ref TNM(TMat_NM.data(), {NAEA,NMO});
      SlaterDeterminantOperations::base::MixedDensityMatrix_noHerm<T>(A,B,std::forward<MatC>(C),res,TNN,TNM,IWORK,WORK,compact);
    }

    template<class integer, class MatA, class MatB, class MatC, class MatQ>
    void MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C, T* res, 
                                    integer* ref, MatQ&& QQ0, bool compact=false) {
      int Nact = hermA.shape()[0];
      int NEL = B.shape()[1];
      int NMO = B.shape()[0];
      assert(hermA.shape()[1]==B.shape()[0]);
      assert(QQ0.shape()[0]==Nact);
      assert(QQ0.shape()[1]==NEL);
      assert(TMat_NN.num_elements() >= NEL*NEL);
      TMatrix_ref TNN(TMat_NN.data(), {NEL,NEL});
      assert(TMat_NM.num_elements() >= Nact*NEL);
      TMatrix_ref TAB(TMat_NM.data(), {Nact,NEL});
      assert(TMat_MM.num_elements() >= NMO*NEL);
      TMatrix_ref TNM(TMat_MM.data(), {NEL,NMO});
      SlaterDeterminantOperations::base::MixedDensityMatrixForWoodbury<T>(hermA,B,std::forward<MatC>(C),res,std::forward<MatQ>(QQ0),ref,TNN,TAB,TNM,IWORK,WORK,compact);
    }

    template<class integer, class MatA, class MatB, class MatC>
    void MixedDensityMatrixFromConfiguration(const MatA& hermA, const MatB& B, MatC&& C, T* res,
                                    integer* ref, bool compact=false) {
      int Nact = hermA.shape()[0];
      int NEL = B.shape()[1];
      int NMO = B.shape()[0];
      assert(hermA.shape()[1]==B.shape()[0]);
      assert(TMat_NN.num_elements() >= NEL*NEL);
      TMatrix_ref TNN(TMat_NN.data(), {NEL,NEL});
      assert(TMat_NM.num_elements() >= Nact*NEL);
      TMatrix_ref TAB(TMat_NM.data(), {Nact,NEL});
      assert(TMat_MM.num_elements() >= NMO*NEL);
      TMatrix_ref TNM(TMat_MM.data(), {NEL,NMO});
      SlaterDeterminantOperations::base::MixedDensityMatrixFromConfiguration<T>(hermA,B,std::forward<MatC>(C),res,ref,TNN,TAB,TNM,IWORK,WORK,compact);
    }

    template<class MatA>
    void Overlap(const MatA& A, T* res) {
      int NAEA = A.shape()[1];
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      TMatrix_ref TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NAEA);
      TMatrix_ref TNN2(TMat_NM.data(), {NAEA,NAEA});
      SlaterDeterminantOperations::base::Overlap_noHerm<T>(A,A,res,TNN,IWORK,TNN2);
    }

    template<class MatA, class MatB>
    void Overlap(const MatA& hermA, const MatB& B, T* res) {
      int NAEA = hermA.shape()[0];
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      TMatrix_ref TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NAEA);
      TMatrix_ref TNN2(TMat_NM.data(), {NAEA,NAEA});
      SlaterDeterminantOperations::base::Overlap<T>(hermA,B,res,TNN,IWORK,TNN2);
    } 

    template<class MatA, class MatB>
    void Overlap_noHerm(const MatA& A, const MatB& B, T* res) {
      int NAEA = A.shape()[1];
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      TMatrix_ref TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NAEA);
      TMatrix_ref TNN2(TMat_NM.data(), {NAEA,NAEA});
      SlaterDeterminantOperations::base::Overlap_noHerm<T>(A,B,res,TNN,IWORK,TNN2);
    }

    // routines for PHMSD
    template<typename integer, class MatA, class MatB, class MatC>
    void OverlapForWoodbury(const MatA& hermA, const MatB& B, T* res, integer* ref, MatC&& QQ0) {
      int Nact = hermA.shape()[0];
      int NEL = B.shape()[1];
      assert(hermA.shape()[1]==B.shape()[0]);
      assert(QQ0.shape()[0]==Nact);  
      assert(QQ0.shape()[1]==NEL);  
      assert(TMat_NN.num_elements() >= NEL*NEL);
      assert(TMat_MM.num_elements() >= Nact*NEL);
      TMatrix_ref TNN(TMat_NN.data(), {NEL,NEL});
      TMatrix_ref TMN(TMat_MM.data(), {Nact,NEL});
      SlaterDeterminantOperations::base::OverlapForWoodbury<T>(hermA,B,res,std::forward<MatC>(QQ0),ref,TNN,TMN,IWORK,WORK);
    }    

    template<class Mat, class MatP1, class MatV>
    void Propagate(Mat&& A, const MatP1& P1, const MatV& V, int order=6) {
      int NMO = A.shape()[0];
      int NAEA = A.shape()[1];
      if(TMat_MN.num_elements() < NMO*NAEA)
        TMat_MN.reextent({NMO,NAEA});
      if(TMat_NM.num_elements() < NMO*NAEA)
        TMat_NM.reextent({NAEA,NMO});
      TMatrix_ref TMN(TMat_MN.data(), {NMO,NAEA});
      TMatrix_ref T1(TMat_NM.data(), {NMO,NAEA});
      TMatrix_ref T2(TMat_MM.data(), {NMO,NAEA});
      ma::product(P1,std::forward<Mat>(A),TMN);
      SlaterDeterminantOperations::base::apply_expM(V,TMN,T1,T2,order);
      ma::product(P1,TMN,std::forward<Mat>(A));
    }

    // need to check if this is equivalent to QR!!!
    template<class Mat>
    void Orthogonalize(Mat&& A, T* res=nullptr) {
      ma::gelqf(std::forward<Mat>(A),TAU,WORK);
      if(res != nullptr) {
        *res = T(1.0); 
        for (int i = 0; i < A.shape()[1]; i++) { 
          if (real(A[i][i]) < 0) 
            IWORK[i]=-1; 
          else 
            IWORK[i]=1; 
          *res *= T(IWORK[i])*A[i][i];
        }
      }
      ma::glq(std::forward<Mat>(A),TAU,WORK);
      for(int i=0; i<A.shape()[0]; ++i)
        for(int j=0; j<A.shape()[1]; ++j)
          A[i][j] *= T(IWORK[j]);
    }

  protected:

    Alloc allocator_;
    IAlloc iallocator_;

    TVector WORK;
    IVector IWORK;

    // Vector used in QR routines 
    TVector TAU;

    // TMat_AB: Local temporary Matrix of dimension [AxB]
    // N: NAEA
    // M: NMO
    TMatrix TMat_NN;
    TMatrix TMat_NM;
    TMatrix TMat_MN;
    TMatrix TMat_MM;
    TMatrix TMat_MM2;
    TMatrix TMat_MM3;

};

}

}

#endif
