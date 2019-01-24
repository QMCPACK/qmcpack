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
    using Base::TMat_NN;
    using Base::TMat_NM;
    using Base::TMat_MN;
    using Base::TMat_MM;
    using Base::TMat_MM2;
    using Base::TMat_MM3;

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
    void MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, T* res, communicator& comm, bool compact=false) {
      int NMO = hermA.size(1);
      int NAEA = hermA.size(0);
      set_shm_buffer(comm,NAEA*(NAEA+NMO));
      assert(SM_TMats->num_elements() >= NAEA*(NAEA+NMO));
      boost::multi::array_ref<T,2> TNN(to_address(SM_TMats->origin()), {NAEA,NAEA});
      boost::multi::array_ref<T,2> TNM(to_address(SM_TMats->origin())+NAEA*NAEA, {NAEA,NMO});  
      SlaterDeterminantOperations::shm::MixedDensityMatrix<T>(hermA,B,std::forward<MatC>(C),res,TNN,TNM,IWORK,WORK,comm,compact);
    }

    template<class integer, class MatA, class MatB, class MatC, class MatQ>
    void MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C, T* res,
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
      SlaterDeterminantOperations::shm::MixedDensityMatrixForWoodbury<T>(hermA,B,std::forward<MatC>(C),res,std::forward<MatQ>(QQ0),ref,TNN,TAB,TNM,IWORK,WORK,comm,compact);
    }

    template<class MatA, class MatB>
    void Overlap(const MatA& hermA, const MatB& B, T* res, communicator& comm) { 
      int NAEA = hermA.size(0);
      set_shm_buffer(comm,2*NAEA*NAEA);
      assert(SM_TMats->num_elements() >= 2*NAEA*NAEA);
      boost::multi::array_ref<T,2> TNN(to_address(SM_TMats->origin()), {NAEA,NAEA});
      boost::multi::array_ref<T,2> TNN2(to_address(SM_TMats->origin())+NAEA*NAEA, {NAEA,NAEA});
      SlaterDeterminantOperations::shm::Overlap<T>(hermA,B,res,TNN,IWORK,TNN2,comm);
    }

    template<typename integer, class MatA, class MatB, class MatC>
    void OverlapForWoodbury(const MatA& hermA, const MatB& B, T* res, integer* ref, MatC&& QQ0, communicator& comm) {
      int Nact = hermA.size(0);
      int NEL = B.size(1);
      assert(hermA.size(1)==B.size(0));
      assert(QQ0.size(0)==Nact);
      assert(QQ0.size(1)==NEL);
      assert(TMat_NN.num_elements() >= NEL*NEL);
      assert(TMat_MM.num_elements() >= Nact*NEL);
      set_shm_buffer(comm,NEL*(Nact+NEL));
      assert(SM_TMats->num_elements() >= NEL*(Nact+NEL));
      boost::multi::array_ref<T,2> TNN(to_address(to_address(SM_TMats->origin())), {NEL,NEL});
      boost::multi::array_ref<T,2> TMN(to_address((to_address(SM_TMats->origin())+NEL*NEL)), {Nact,NEL});
      SlaterDeterminantOperations::shm::OverlapForWoodbury<T>(hermA,B,res,std::forward<MatC>(QQ0),ref,TNN,TMN,IWORK,WORK,comm);
    }

    template<class Mat, class MatP1, class MatV>
    void Propagate(Mat&& A, const MatP1& P1, const MatV& V, communicator& comm, int order=6) {
      int NMO = A.size(0);
      int NAEA = A.size(1);
      set_shm_buffer(comm,3*NAEA*NMO);
      assert(SM_TMats->num_elements() >= 3*NAEA*NAEA);
      boost::multi::array_ref<T,2> T0(to_address(SM_TMats->origin()), {NMO,NAEA});
      boost::multi::array_ref<T,2> T1(to_address(SM_TMats->origin())+NMO*NAEA, {NMO,NAEA});
      boost::multi::array_ref<T,2> T2(to_address(SM_TMats->origin())+2*NMO*NAEA, {NMO,NAEA});
      if(comm.root()) 
        ma::product(P1,std::forward<Mat>(A),T0);
      comm.barrier();
      SlaterDeterminantOperations::shm::apply_expM(V,T0,T1,T2,comm,order);
      comm.barrier();
      if(comm.root()) 
        ma::product(P1,T0,std::forward<Mat>(A));
      comm.barrier();
    }

  private:

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
