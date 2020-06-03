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
#include "AFQMC/Utilities/Utils.hpp"
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

template<typename Type, class buffer_generator_type> 
class SlaterDetOperations_base 
{
  public:

    using T = Type;
    using R = typename remove_complex<T>::value_type; 

    using buffer_type_I = typename buffer_generator_type::template allocator<int>;
    using buffer_type_R = typename buffer_generator_type::template allocator<R>;
    using buffer_type_T = typename buffer_generator_type::template allocator<T>;

    using IVector = boost::multi::static_array<int,1,buffer_type_I>; 
    using RVector = boost::multi::static_array<R,1,buffer_type_R>; 
    using TVector = boost::multi::static_array<T,1,buffer_type_T>; 
    using TMatrix = boost::multi::static_array<T,2,buffer_type_T>; 

    SlaterDetOperations_base(buffer_generator_type* bgen):
      buffer_generator(bgen)
    {
    }

    SlaterDetOperations_base(int NMO, int NAEA, buffer_generator_type* bgen):
      buffer_generator(bgen)
    {

      // only used to determine size of WORK  
      TMatrix TNN({NAEA,NAEA},buffer_generator->template get_allocator<T>());
      TMatrix TNM({NAEA,NMO},buffer_generator->template get_allocator<T>());
      TMatrix TMN({NMO,NAEA},buffer_generator->template get_allocator<T>());

      // reserve enough space in lapack's work array
      // Make sure it is large enough for:
      // 1. getri( TNN )
      //  2. geqrf( TNM )
      work_size = std::max(ma::getri_optimal_workspace_size(TNN),
                                  ma::geqrf_optimal_workspace_size(TNM));  
      //  3. gqr( TNM )
      work_size = std::max(work_size, ma::gqr_optimal_workspace_size(TNM) );
      //  4. gelqf( TMN )
      work_size = std::max(work_size, ma::gelqf_optimal_workspace_size(TMN) ); 
      //  5. glq( TMN )
      work_size = std::max(work_size, ma::glq_optimal_workspace_size(TMN) ); 
      //  6. trf( TNN )
      work_size = std::max(work_size, ma::getrf_optimal_workspace_size(TNN) ); 
      //  7. svd( TNN )
      work_size = std::max(work_size, ma::gesvd_optimal_workspace_size(TNN) ); 
      work_size = std::max(work_size, NMO); 
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
      TMatrix TNN({NAEA,NAEA},buffer_generator->template get_allocator<T>());
      TMatrix TNM({NAEA,NMO},buffer_generator->template get_allocator<T>());
      TVector WORK(iextensions<1u>{work_size},buffer_generator->template get_allocator<T>());
      IVector IWORK(iextensions<1u>{NMO+1},buffer_generator->template get_allocator<int>());
      return SlaterDeterminantOperations::base::MixedDensityMatrix<T>(hermA,B,std::forward<MatC>(C),LogOverlapFactor,TNN,TNM,IWORK,WORK,compact,herm);
    }

    template<class MatA, class MatC>
    T MixedDensityMatrix(const MatA& A, MatC&& C, T LogOverlapFactor, bool compact=false) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      TMatrix TNN({NAEA,NAEA},buffer_generator->template get_allocator<T>());
      TMatrix TNM({NAEA,NMO},buffer_generator->template get_allocator<T>());
      TVector WORK(iextensions<1u>{work_size},buffer_generator->template get_allocator<T>());
      IVector IWORK(iextensions<1u>{NMO+1},buffer_generator->template get_allocator<int>());
      return SlaterDeterminantOperations::base::MixedDensityMatrix<T>(A,A,std::forward<MatC>(C),LogOverlapFactor,TNN,TNM,IWORK,WORK,compact,false);
    }

    template<class MatA, class MatB, class MatC>
    T MixedDensityMatrix_noHerm(const MatA& A, const MatB& B, MatC&& C, T LogOverlapFactor, bool compact=false, bool useSVD = false) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      if(useSVD) {
        TMatrix TNN1({NAEA,NAEA},buffer_generator->template get_allocator<T>());
        TMatrix TNN2({NAEA,NAEA},buffer_generator->template get_allocator<T>());
        TMatrix TNM({NAEA,NMO},buffer_generator->template get_allocator<T>());
        TMatrix TMN({NMO,NAEA},buffer_generator->template get_allocator<T>());
        TVector WORK(iextensions<1u>{work_size},buffer_generator->template get_allocator<T>());
        RVector RWORK(iextensions<1u>{6*NAEA+1},buffer_generator->template get_allocator<R>());
        IVector IWORK(iextensions<1u>{NMO+1},buffer_generator->template get_allocator<int>());
        return SlaterDeterminantOperations::base::MixedDensityMatrix_noHerm_wSVD<T>(A,B,std::forward<MatC>(C),LogOverlapFactor,RWORK,TNN1,TNN2,TMN,TNM,IWORK,WORK,compact);
      } else {
        TMatrix TNN({NAEA,NAEA},buffer_generator->template get_allocator<T>());
        TMatrix TNM({NAEA,NMO},buffer_generator->template get_allocator<T>());
        TVector WORK(iextensions<1u>{work_size},buffer_generator->template get_allocator<T>());
        IVector IWORK(iextensions<1u>{NMO+1},buffer_generator->template get_allocator<int>());
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
      TMatrix TNN({NEL,NEL},buffer_generator->template get_allocator<T>());
      TMatrix TAB({Nact,NEL},buffer_generator->template get_allocator<T>());
      TMatrix TNM({NEL,NMO},buffer_generator->template get_allocator<T>());
      TVector WORK(iextensions<1u>{work_size},buffer_generator->template get_allocator<T>());
      IVector IWORK(iextensions<1u>{NMO+1},buffer_generator->template get_allocator<int>());
      return SlaterDeterminantOperations::base::MixedDensityMatrixForWoodbury<T>(hermA,B,std::forward<MatC>(C),LogOverlapFactor,std::forward<MatQ>(QQ0),ref,TNN,TAB,TNM,IWORK,WORK,compact);
    }

    template<class integer, class MatA, class MatB, class MatC>
    T MixedDensityMatrixFromConfiguration(const MatA& hermA, const MatB& B, MatC&& C, T LogOverlapFactor, 
                                    integer* ref, bool compact=false) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      int NMO = B.size(0);
      assert(hermA.size(1)==B.size(0));
      TMatrix TNN({NEL,NEL},buffer_generator->template get_allocator<T>());
      TMatrix TAB({Nact,NEL},buffer_generator->template get_allocator<T>());
      TMatrix TNM({NEL,NMO},buffer_generator->template get_allocator<T>());
      TVector WORK(iextensions<1u>{work_size},buffer_generator->template get_allocator<T>());
      IVector IWORK(iextensions<1u>{NMO+1},buffer_generator->template get_allocator<int>());
      return SlaterDeterminantOperations::base::MixedDensityMatrixFromConfiguration<T>(hermA,B,std::forward<MatC>(C),LogOverlapFactor,ref,TNN,TAB,TNM,IWORK,WORK,compact);
    }

    template<class MatA>
    T Overlap(const MatA& A, T LogOverlapFactor) {
      int NAEA = A.size(1);
      TMatrix TNN({NAEA,NAEA},buffer_generator->template get_allocator<T>());
      TMatrix TNN2({NAEA,NAEA},buffer_generator->template get_allocator<T>());
      IVector IWORK(iextensions<1u>{NAEA+1},buffer_generator->template get_allocator<int>());
      return SlaterDeterminantOperations::base::Overlap<T>(A,A,LogOverlapFactor,TNN,IWORK,TNN2,false);
    }

    template<class MatA, class MatB>
    T Overlap(const MatA& hermA, const MatB& B, T LogOverlapFactor, bool herm=true) {
      int NAEA = (herm?hermA.size(0):hermA.size(1));
      TMatrix TNN({NAEA,NAEA},buffer_generator->template get_allocator<T>());
      TMatrix TNN2({NAEA,NAEA},buffer_generator->template get_allocator<T>());
      IVector IWORK(iextensions<1u>{NAEA+1},buffer_generator->template get_allocator<int>());
      return SlaterDeterminantOperations::base::Overlap<T>(hermA,B,LogOverlapFactor,TNN,IWORK,TNN2,herm);
    } 

    template<class MatA, class MatB>
    T Overlap_noHerm(const MatA& A, const MatB& B, T LogOverlapFactor) {
      int NAEA = A.size(1);
      TMatrix TNN({NAEA,NAEA},buffer_generator->template get_allocator<T>());
      TMatrix TNN2({NAEA,NAEA},buffer_generator->template get_allocator<T>());
      IVector IWORK(iextensions<1u>{NAEA+1},buffer_generator->template get_allocator<int>());
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
      TMatrix TNN({NEL,NEL},buffer_generator->template get_allocator<T>());
      TMatrix TMN({Nact,NEL},buffer_generator->template get_allocator<T>());
      TVector WORK(iextensions<1u>{work_size},buffer_generator->template get_allocator<T>());
      IVector IWORK(iextensions<1u>{Nact+1},buffer_generator->template get_allocator<int>());
      return SlaterDeterminantOperations::base::OverlapForWoodbury<T>(hermA,B,LogOverlapFactor,std::forward<MatC>(QQ0),ref,TNN,TMN,IWORK,WORK);
    }    

    template<class Mat, class MatP1, class MatV>
    void Propagate(Mat&& A, const MatP1& P1, const MatV& V, int order=6, char TA='N') {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      TMatrix TMN({NMO,NAEA},buffer_generator->template get_allocator<T>());
      TMatrix T1({NMO,NAEA},buffer_generator->template get_allocator<T>());
      TMatrix T2({NMO,NAEA},buffer_generator->template get_allocator<T>());
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
      TMatrix AT({NAEA,NMO},buffer_generator->template get_allocator<T>());
      TVector scl(iextensions<1u>{NMO},buffer_generator->template get_allocator<T>());
      TVector TAU(iextensions<1u>{NMO},buffer_generator->template get_allocator<T>());
      TVector WORK(iextensions<1u>{work_size},buffer_generator->template get_allocator<T>());
      IVector IWORK(iextensions<1u>{NMO+1},buffer_generator->template get_allocator<int>());
      ma::transpose(A,AT);   
      ma::geqrf(AT,TAU,WORK);
      using ma::determinant_from_geqrf;
      using ma::scale_columns;  
      T res = determinant_from_geqrf(AT.size(0),AT.origin(),AT.stride(0),scl.origin(),LogOverlapFactor);
      ma::gqr(AT,TAU,WORK);
      ma::transpose(AT,A);   
      scale_columns(A.size(0),A.size(1),A.origin(),A.stride(0),scl.origin());
#else
      int NMO = A.size(0);
      TVector TAU(iextensions<1u>{NMO},buffer_generator->template get_allocator<T>());
      TVector WORK(iextensions<1u>{work_size},buffer_generator->template get_allocator<T>());
      IVector IWORK(iextensions<1u>{NMO+1},buffer_generator->template get_allocator<int>());
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

    buffer_generator_type* buffer_generator;
    //host_allocator_generator_type* buffer_generator;

    // keeping this to avoid having to know which routines are called in the lower level
    //  and let's me use static arrays
    int work_size=0;

};

}

}

#endif
