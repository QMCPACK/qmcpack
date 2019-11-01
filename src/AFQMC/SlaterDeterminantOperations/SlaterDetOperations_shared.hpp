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

#ifndef QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_SHARED_HPP
#define QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_SHARED_HPP

#include<fstream>

#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/csr_blas.hpp"
#include "AFQMC/SlaterDeterminantOperations/mixed_density_matrix.hpp"
#include "AFQMC/SlaterDeterminantOperations/apply_expM.hpp"

#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations_base.hpp"

#include "mpi3/shared_communicator.hpp"
#include "type_traits/scalar_traits.h"

namespace qmcplusplus
{

namespace afqmc
{

template<class T = ComplexType>
class SlaterDetOperations_shared : public SlaterDetOperations_base<std::allocator<T>>
{
  public:

    using Base = SlaterDetOperations_base<std::allocator<T>>;
    using Alloc = std::allocator<T>; 
    using pointer = typename Base::pointer;
    using const_pointer = typename Base::const_pointer;
    using IAlloc = typename Base::IAlloc;
    using communicator = boost::mpi3::shared_communicator;
    using shmTVector = boost::multi::array<T,1,shared_allocator<T>>;  

    using IVector = typename Base::IVector;
    using TVector = typename Base::TVector;
    using TMatrix = typename Base::TMatrix;

    using Base::MixedDensityMatrix;
    using Base::MixedDensityMatrixForWoodbury;
    using Base::MixedDensityMatrix_noHerm;
    using Base::MixedDensityMatrixFromConfiguration;
    using Base::Overlap;
    using Base::Overlap_noHerm;
    using Base::OverlapForWoodbury;
    using Base::Propagate;
    using Base::Orthogonalize;
    using Base::IWORK;
    using Base::WORK;

    SlaterDetOperations_shared():
      SlaterDetOperations_base<Alloc>(Alloc{}),
      SM_TMats(nullptr)
    {
    }

    SlaterDetOperations_shared(int NMO, int NAEA):
      SlaterDetOperations_base<Alloc>(NMO,NAEA,Alloc{}),
      SM_TMats(nullptr)
    {
    }

    ~SlaterDetOperations_shared() {}

    SlaterDetOperations_shared(const SlaterDetOperations_shared& other) = delete;
    SlaterDetOperations_shared(SlaterDetOperations_shared&& other) = default;
    SlaterDetOperations_shared& operator=(const SlaterDetOperations_shared& other) = delete;
    SlaterDetOperations_shared& operator=(SlaterDetOperations_shared&& other) = default;

    // C must live in shared memory for this routine to work as expected
    template<class MatA, class MatB, class MatC>
    T MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, T LogOverlapFactor, communicator& comm, bool compact=false, bool herm=true) {
      int NMO = (herm?hermA.size(1):hermA.size(0));
      int NAEA = (herm?hermA.size(0):hermA.size(1));
      set_shm_buffer(comm,NAEA*(NAEA+NMO));
      assert(SM_TMats->num_elements() >= NAEA*(NAEA+NMO));
      boost::multi::array_ref<T,2> TNN(to_address(SM_TMats->origin()), {NAEA,NAEA});
      boost::multi::array_ref<T,2> TNM(to_address(SM_TMats->origin())+NAEA*NAEA, {NAEA,NMO});  
      return SlaterDeterminantOperations::shm::MixedDensityMatrix<T>(hermA,B,std::forward<MatC>(C),LogOverlapFactor,TNN,TNM,IWORK,WORK,comm,compact,herm);
    }

    template<class integer, class MatA, class MatB, class MatC, class MatQ>
    T MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C, T LogOverlapFactor, 
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
      boost::multi::array_ref<T,2> TNN(to_address(SM_TMats->origin()), {NEL,NEL});
      cnt+=TNN.num_elements();
      boost::multi::array_ref<T,2> TAB(to_address(SM_TMats->origin())+cnt, {Nact,NEL});
      cnt+=TAB.num_elements();
      boost::multi::array_ref<T,2> TNM(to_address(SM_TMats->origin())+cnt, {NEL,NMO});
      return SlaterDeterminantOperations::shm::MixedDensityMatrixForWoodbury<T>(hermA,B,std::forward<MatC>(C),LogOverlapFactor,std::forward<MatQ>(QQ0),ref,TNN,TAB,TNM,IWORK,WORK,comm,compact);
    }

    template<class MatA, class MatB>
    T Overlap(const MatA& hermA, const MatB& B, T LogOverlapFactor, communicator& comm, bool herm=true) { 
      int NAEA = (herm?hermA.size(0):hermA.size(1));
      set_shm_buffer(comm,2*NAEA*NAEA);
      assert(SM_TMats->num_elements() >= 2*NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN(to_address(SM_TMats->origin()), {NAEA,NAEA});
      boost::multi::array_ref<T,2> TNN2(to_address(SM_TMats->origin())+NAEA*NAEA, {NAEA,NAEA});
      return SlaterDeterminantOperations::shm::Overlap<T>(hermA,B,LogOverlapFactor,TNN,IWORK,TNN2,comm,herm);
    }

    template<typename integer, class MatA, class MatB, class MatC>
    T OverlapForWoodbury(const MatA& hermA, const MatB& B, T LogOverlapFactor, integer* ref, MatC&& QQ0, communicator& comm) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      assert(hermA.size(1)==B.size(0));
      assert(QQ0.size(0)==Nact);
      assert(QQ0.size(1)==NEL);
      set_shm_buffer(comm,NEL*(Nact+NEL));
      assert(SM_TMats->num_elements() >= NEL*(Nact+NEL));
      boost::multi::array_ref<T,2> TNN(to_address(SM_TMats->origin()), {NEL,NEL});
      boost::multi::array_ref<T,2> TMN(to_address(SM_TMats->origin())+NEL*NEL, {Nact,NEL});
      return SlaterDeterminantOperations::shm::OverlapForWoodbury<T>(hermA,B,LogOverlapFactor,std::forward<MatC>(QQ0),ref,TNN,TMN,IWORK,WORK,comm);
    }

    template<class Mat, class MatP1, class MatV>
    void Propagate(Mat&& A, const MatP1& P1, const MatV& V, communicator& comm, int order=6, char TA='N') {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      set_shm_buffer(comm,3*NAEA*NMO);
      assert(SM_TMats->num_elements() >= 3*NAEA*NAEA);
      boost::multi::array_ref<T,2> T0(to_address(SM_TMats->origin()), {NMO,NAEA});
      boost::multi::array_ref<T,2> T1(to_address(SM_TMats->origin())+NMO*NAEA, {NMO,NAEA});
      boost::multi::array_ref<T,2> T2(to_address(SM_TMats->origin())+2*NMO*NAEA, {NMO,NAEA});
      using ma::T;
      using ma::H;
      if(comm.root()) { 
        if(TA=='H' || TA=='h') 
          ma::product(ma::H(P1),std::forward<Mat>(A),T0);
        else if(TA=='T' || TA=='t') 
          ma::product(ma::T(P1),std::forward<Mat>(A),T0);
        else 
          ma::product(P1,std::forward<Mat>(A),T0);
      }  
      comm.barrier();
      SlaterDeterminantOperations::shm::apply_expM(V,T0,T1,T2,comm,order,TA);
      comm.barrier();
      if(comm.root()) {
        if(TA=='H' || TA=='h') 
          ma::product(ma::H(P1),std::forward<Mat>(A),T0);
        else if(TA=='T' || TA=='t') 
          ma::product(ma::T(P1),std::forward<Mat>(A),T0);
        else 
          ma::product(P1,std::forward<Mat>(A),T0);
      }
      comm.barrier();
    }

    // C[nwalk, M, N]
    template<class MatA, class MatB, class MatC, class TVec>
    void BatchedMixedDensityMatrix(const MatA& hermA, std::vector<MatB> &Bi, MatC&& C, T LogOverlapFactor, TVec&& ovlp, bool compact=false, bool herm=true) {
      APP_ABORT(" Error: Batched routines not compatible with SlaterDetOperations_shared::BatchedMixedDensityMatrix \n");
    }

    template<class MatA, class MatB, class MatC, class TVec>
    void BatchedDensityMatrices(const std::vector<MatA>& hermA, std::vector<MatB> &Bi, std::vector<MatC>& C, T LogOverlapFactor, TVec&& ovlp, bool compact=false, bool herm=true) {
      APP_ABORT(" Error: Batched routines not compatible with SlaterDetOperations_shared::BatchedMixedDensityMatrix \n");
    }

    template<class MatA, class MatB, class TVec>
    void BatchedOverlap(const MatA& hermA, std::vector<MatB> &Bi, T LogOverlapFactor, TVec&& ovlp, bool herm=true) {
      APP_ABORT(" Error: Batched routines not compatible with SlaterDetOperations_shared::BatchedOverlap \n");
    }

    template<class MatA, class MatP1, class MatV>
    void BatchedPropagate(std::vector<MatA> &Ai, const MatP1& P1, const MatV& V, int order=6, char TA='N') {
      APP_ABORT(" Error: Batched routines not compatible with SlaterDetOperations_shared::BatchedPropagate \n");
    }

    template<class MatA>
    void BatchedOrthogonalize(std::vector<MatA> &Ai, T LogOverlapFactor) {
      APP_ABORT(" Error: Batched routines not compatible with SlaterDetOperations_shared::BatchedOrthogonalize \n");
    }

    template<class MatA, class PTR>
    void BatchedOrthogonalize(std::vector<MatA> &Ai, T LogOverlapFactor, PTR detR) {
      APP_ABORT(" Error: Batched routines not compatible with SlaterDetOperations_shared::BatchedOrthogonalize \n");
    }

  protected:

    // shm temporary matrices
    std::unique_ptr<shmTVector> SM_TMats;

    void set_shm_buffer(communicator& comm, size_t N) {
      if( SM_TMats == nullptr || SM_TMats->get_allocator() != shared_allocator<T>{comm} ) { 
        SM_TMats = std::move(std::make_unique<shmTVector>(iextensions<1u>{N},shared_allocator<T>{comm}));
      } else if(SM_TMats->num_elements() < N) 
        SM_TMats = std::move(std::make_unique<shmTVector>(iextensions<1u>{N},shared_allocator<T>{comm}));
    }

};

}

}

#endif
