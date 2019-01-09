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

#ifndef QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_HPP
#define QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_HPP

#include<fstream>

#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/csr_blas.hpp"
#include "AFQMC/SlaterDeterminantOperations/mixed_density_matrix.hpp"
#include "AFQMC/SlaterDeterminantOperations/apply_expM.hpp"

#include "mpi3/shared_communicator.hpp"
#include "type_traits/scalar_traits.h"

namespace qmcplusplus
{

namespace afqmc
{

template<class T = ComplexType> 
class SlaterDetOperations 
{
  using communicator = boost::mpi3::shared_communicator;
  using IVector = boost::multi::array<int,1>;  
  using TVector = boost::multi::array<T,1>;  
  using TMatrix = boost::multi::array<T,2>;  
  using shmTVector = boost::multi::array<T,1,shared_allocator<T>>;  

  public:

    SlaterDetOperations(int NMO, int NAEA):
        SM_TMats(nullptr)
    {
      // IWORK: integer buffer for getri/getrf, which expects NMO+1 
      IWORK.reextent( extensions<1u>{NMO+1});

      // local temporary storage
      TMat_NM.reextent({NAEA,NMO});
      TMat_MN.reextent({NMO,NAEA});
      TMat_NN.reextent({NAEA,NAEA});
      TMat_MM.reextent({NMO,NMO});
      TMat_MM2.reextent({NMO,NMO});
      TMat_MM3.reextent({NMO,NMO});

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

      // TAU: T used in QR routines 
      TAU.reextent(extensions<1u>{NMO});
    }

    ~SlaterDetOperations() {}

    SlaterDetOperations(const SlaterDetOperations& other) = delete;
    SlaterDetOperations(SlaterDetOperations&& other) = default;
    SlaterDetOperations& operator=(const SlaterDetOperations& other) = delete;
    SlaterDetOperations& operator=(SlaterDetOperations&& other) = default;

    template<class MatA, class MatB, class MatC>
    T MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, bool compact=false) {
      int NMO = hermA.size(1);
      int NAEA = hermA.size(0);
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NMO);
      boost::multi::array_ref<T,2> TNM(TMat_NM.data(), {NAEA,NMO});
      return SlaterDeterminantOperations::base::MixedDensityMatrix<T>(hermA,B,std::forward<MatC>(C),TNN,TNM,IWORK,WORK,compact);
    }

    template<class MatA, class MatC>
    T MixedDensityMatrix(const MatA& A, MatC&& C, bool compact=false) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NMO);
      boost::multi::array_ref<T,2> TNM(TMat_NM.data(), {NAEA,NMO});
      return SlaterDeterminantOperations::base::MixedDensityMatrix_noHerm<T>(A,A,std::forward<MatC>(C),TNN,TNM,IWORK,WORK,compact);
    }

    template<class MatA, class MatB, class MatC>
    T MixedDensityMatrix_noHerm(const MatA& A, const MatB& B, MatC&& C, bool compact=false) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NMO);
      boost::multi::array_ref<T,2> TNM(TMat_NM.data(), {NAEA,NMO});
      return SlaterDeterminantOperations::base::MixedDensityMatrix_noHerm<T>(A,B,std::forward<MatC>(C),TNN,TNM,IWORK,WORK,compact);
    }

    // C must live in shared memory for this routine to work as expected
    template<class MatA, class MatB, class MatC>
    T MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, communicator& comm, bool compact=false) {
      int NMO = hermA.size(1);
      int NAEA = hermA.size(0);
      set_shm_buffer(comm,NAEA*(NAEA+NMO));
      assert(SM_TMats->num_elements() >= NAEA*(NAEA+NMO));
      boost::multi::array_ref<T,2> TNN(std::addressof(*SM_TMats->origin()), {NAEA,NAEA});
      boost::multi::array_ref<T,2> TNM(std::addressof(*SM_TMats->origin())+NAEA*NAEA, {NAEA,NMO});  
      return SlaterDeterminantOperations::shm::MixedDensityMatrix<T>(hermA,B,std::forward<MatC>(C),TNN,TNM,IWORK,WORK,comm,compact);
    }

    template<class integer, class MatA, class MatB, class MatC, class MatQ>
    T MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C, 
                                    integer* ref, MatQ&& QQ0, bool compact=false) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      int NMO = B.size(0);
      assert(hermA.size(1)==B.size(0));
      assert(QQ0.size(0)==Nact);
      assert(QQ0.size(1)==NEL);
      assert(TMat_NN.num_elements() >= NEL*NEL);
      boost::multi::array_ref<T,2> TNN(TMat_NN.data(), {NEL,NEL});
      assert(TMat_NM.num_elements() >= Nact*NEL);
      boost::multi::array_ref<T,2> TAB(TMat_NM.data(), {Nact,NEL});
      assert(TMat_MM.num_elements() >= NMO*NEL);
      boost::multi::array_ref<T,2> TNM(TMat_MM.data(), {NEL,NMO});
      return SlaterDeterminantOperations::base::MixedDensityMatrixForWoodbury<T>(hermA,B,std::forward<MatC>(C),std::forward<MatQ>(QQ0),ref,TNN,TAB,TNM,IWORK,WORK,compact);
    }

    template<class integer, class MatA, class MatB, class MatC, class MatQ>
    T MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C,
                                    integer* ref, MatQ&& QQ0, communicator& comm, bool compact=false) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      int NMO = B.size(0);
      assert(hermA.size(1)==B.size(0));
      assert(QQ0.size(0)==Nact);
      assert(QQ0.size(1)==NEL);

      set_shm_buffer(comm,NEL*(NEL+Nact+NMO));
      assert(SM_TMats->num_elements() >= NEL*(NEL+Nact+NMO));
      size_t cnt=0;
      boost::multi::array_ref<T,2> TNN(std::addressof(*SM_TMats->origin()), {NEL,NEL});
      cnt+=TNN.num_elements();
      boost::multi::array_ref<T,2> TAB(std::addressof(*SM_TMats->origin())+cnt, {Nact,NEL});
      cnt+=TAB.num_elements();
      boost::multi::array_ref<T,2> TNM(std::addressof(*SM_TMats->origin())+cnt, {NEL,NMO});
      return SlaterDeterminantOperations::shm::MixedDensityMatrixForWoodbury<T>(hermA,B,std::forward<MatC>(C),std::forward<MatQ>(QQ0),ref,TNN,TAB,TNM,IWORK,WORK,comm,compact);
    }

    template<class integer, class MatA, class MatB, class MatC>
    T MixedDensityMatrixFromConfiguration(const MatA& hermA, const MatB& B, MatC&& C,
                                    integer* ref, bool compact=false) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      int NMO = B.size(0);
      assert(hermA.size(1)==B.size(0));
      assert(TMat_NN.num_elements() >= NEL*NEL);
      boost::multi::array_ref<T,2> TNN(TMat_NN.data(), {NEL,NEL});
      assert(TMat_NM.num_elements() >= Nact*NEL);
      boost::multi::array_ref<T,2> TAB(TMat_NM.data(), {Nact,NEL});
      assert(TMat_MM.num_elements() >= NMO*NEL);
      boost::multi::array_ref<T,2> TNM(TMat_MM.data(), {NEL,NMO});
      return SlaterDeterminantOperations::base::MixedDensityMatrixFromConfiguration<T>(hermA,B,std::forward<MatC>(C),ref,TNN,TAB,TNM,IWORK,WORK,compact);
    }

    template<class MatA>
    T Overlap(const MatA& A) {
      int NAEA = A.size(1);
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN2(TMat_NM.data(), {NAEA,NAEA});
      return SlaterDeterminantOperations::base::Overlap_noHerm<T>(A,A,TNN,IWORK,TNN2);
    }

    template<class MatA, class MatB>
    T Overlap(const MatA& hermA, const MatB& B) {
      int NAEA = hermA.size(0);
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN2(TMat_NM.data(), {NAEA,NAEA});
      return SlaterDeterminantOperations::base::Overlap<T>(hermA,B,TNN,IWORK,TNN2);
    } 

    template<class MatA, class MatB>
    T Overlap_noHerm(const MatA& A, const MatB& B) {
      int NAEA = A.size(1);
      assert(TMat_NN.num_elements() >= NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN(TMat_NN.data(), {NAEA,NAEA});
      assert(TMat_NM.num_elements() >= NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN2(TMat_NM.data(), {NAEA,NAEA});
      return SlaterDeterminantOperations::base::Overlap_noHerm<T>(A,B,TNN,IWORK,TNN2);
    }

    template<class MatA, class MatB>
    T Overlap(const MatA& hermA, const MatB& B, communicator& comm) { 
      int NAEA = hermA.size(0);
      set_shm_buffer(comm,2*NAEA*NAEA);
      assert(SM_TMats->num_elements() >= 2*NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN(std::addressof(*SM_TMats->origin()), {NAEA,NAEA});
      boost::multi::array_ref<T,2> TNN2(std::addressof(*SM_TMats->origin())+NAEA*NAEA, {NAEA,NAEA});
      return SlaterDeterminantOperations::shm::Overlap<T>(hermA,B,TNN,IWORK,TNN2,comm);
    }

    // routines for PHMSD
    template<typename integer, class MatA, class MatB, class MatC>
    T OverlapForWoodbury(const MatA& hermA, const MatB& B, integer* ref, MatC&& QQ0) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      assert(hermA.size(1)==B.size(0));
      assert(QQ0.size(0)==Nact);  
      assert(QQ0.size(1)==NEL);  
      assert(TMat_NN.num_elements() >= NEL*NEL);
      assert(TMat_MM.num_elements() >= Nact*NEL);
      boost::multi::array_ref<T,2> TNN(TMat_NN.data(), {NEL,NEL});
      boost::multi::array_ref<T,2> TMN(TMat_MM.data(), {Nact,NEL});
      return SlaterDeterminantOperations::base::OverlapForWoodbury<T>(hermA,B,std::forward<MatC>(QQ0),ref,TNN,TMN,IWORK,WORK);
    }    

    template<typename integer, class MatA, class MatB, class MatC>
    T OverlapForWoodbury(const MatA& hermA, const MatB& B, integer* ref, MatC&& QQ0, communicator& comm) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      assert(hermA.size(1)==B.size(0));
      assert(QQ0.size(0)==Nact);
      assert(QQ0.size(1)==NEL);
      assert(TMat_NN.num_elements() >= NEL*NEL);
      assert(TMat_MM.num_elements() >= Nact*NEL);
      set_shm_buffer(comm,NEL*(Nact+NEL));
      assert(SM_TMats->num_elements() >= NEL*(Nact+NEL));
      boost::multi::array_ref<T,2> TNN(std::addressof(*std::addressof(*SM_TMats->origin())), {NEL,NEL});
      boost::multi::array_ref<T,2> TMN(std::addressof(*(std::addressof(*SM_TMats->origin())+NEL*NEL)), {Nact,NEL});
      return SlaterDeterminantOperations::shm::OverlapForWoodbury<T>(hermA,B,std::forward<MatC>(QQ0),ref,TNN,TMN,IWORK,WORK,comm);
    }

    template<class Mat, class MatP1, class MatV>
    void Propagate(Mat&& A, const MatP1& P1, const MatV& V, int order=6) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      if(TMat_MN.num_elements() < NMO*NAEA)
        TMat_MN.reextent({NMO,NAEA});
      if(TMat_NM.num_elements() < NMO*NAEA)
        TMat_NM.reextent({NAEA,NMO});
      boost::multi::array_ref<T,2> TMN(TMat_MN.data(), {NMO,NAEA});
      boost::multi::array_ref<T,2> T1(TMat_NM.data(), {NMO,NAEA});
      boost::multi::array_ref<T,2> T2(TMat_MM.data(), {NMO,NAEA});
      ma::product(P1,std::forward<Mat>(A),TMN);
      SlaterDeterminantOperations::base::apply_expM(V,TMN,T1,T2,order);
      ma::product(P1,TMN,std::forward<Mat>(A));
    }

    template<class Mat, class MatP1, class MatV>
    void Propagate(Mat&& A, const MatP1& P1, const MatV& V, communicator& comm, int order=6) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      set_shm_buffer(comm,3*NAEA*NMO);
      assert(SM_TMats->num_elements() >= 3*NAEA*NAEA);
      boost::multi::array_ref<T,2> T0(std::addressof(*SM_TMats->origin()), {NMO,NAEA});
      boost::multi::array_ref<T,2> T1(std::addressof(*SM_TMats->origin())+NMO*NAEA, {NMO,NAEA});
      boost::multi::array_ref<T,2> T2(std::addressof(*SM_TMats->origin())+2*NMO*NAEA, {NMO,NAEA});
      if(comm.root()) 
        ma::product(P1,std::forward<Mat>(A),T0);
      comm.barrier();
      SlaterDeterminantOperations::shm::apply_expM(V,T0,T1,T2,comm,order);
      comm.barrier();
      if(comm.root()) 
        ma::product(P1,T0,std::forward<Mat>(A));
      comm.barrier();
    }

    // need to check if this is equivalent to QR!!!
    template<class Mat>
    T Orthogonalize(Mat&& A) {
      T detR = T(1.0);
      ma::gelqf(std::forward<Mat>(A),TAU,WORK);
      for (int i = 0; i < A.size(1); i++) { 
        if (real(A[i][i]) < 0) 
          IWORK[i]=-1; 
        else 
          IWORK[i]=1; 
        detR *= T(IWORK[i])*A[i][i];
      }
      ma::glq(std::forward<Mat>(A),TAU,WORK);
      for(int i=0; i<A.size(0); ++i)
        for(int j=0; j<A.size(1); ++j)
          A[i][j] *= T(IWORK[j]);
      return detR;
    }

  private:

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

    // shm temporary matrices
    std::unique_ptr<shmTVector> SM_TMats;

    // somewhat sensitive piece of code
    void set_shm_buffer(communicator& comm, size_t N) {
      // since there is no way to extract the communicator from SM_TMats  
      if( SM_TMats == nullptr || SM_TMats->get_allocator() != shared_allocator<T>{comm} ) { 
        SM_TMats = std::move(std::make_unique<shmTVector>(extensions<1u>{N},shared_allocator<T>{comm}));
      } else if(SM_TMats->num_elements() < N) 
        SM_TMats->reextent(extensions<1u>{N});
    }
};

}

}

#endif
