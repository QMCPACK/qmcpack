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

#ifndef QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_SERIAL_HPP
#define QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_SERIAL_HPP

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

// Implementation that doesn't use shared memory
// This version is designed for GPU and/or threading with OpenMP
template<class AllocType> 
class SlaterDetOperations_serial : public SlaterDetOperations_base<AllocType> 
{
  public:

    using communicator = boost::mpi3::shared_communicator;
    using Base = SlaterDetOperations_base<AllocType>; 
    using Alloc = typename Base::Alloc;
    using T = typename Base::T;
    using pointer = typename Base::pointer;
    using const_pointer = typename Base::const_pointer;
    using IAlloc = typename Base::IAlloc;

    using IVector = typename Base::IVector;
    using TVector = typename Base::TVector;
    using TMatrix = typename Base::TMatrix;
    using TMatrix_ref = typename Base::TMatrix_ref;
    using TTensor_ref = typename Base::TTensor_ref;

    using Base::MixedDensityMatrix;
    using Base::MixedDensityMatrixForWoodbury;
    using Base::MixedDensityMatrix_noHerm;
    using Base::MixedDensityMatrixFromConfiguration;
    using Base::Overlap;
    using Base::Overlap_noHerm;
    using Base::OverlapForWoodbury;
    using Base::Propagate;
    using Base::Orthogonalize;

    SlaterDetOperations_serial(AllocType alloc_={}):
      SlaterDetOperations_base<AllocType>(alloc_)
    {
    }

    SlaterDetOperations_serial(int NMO, int NAEA, AllocType alloc_={}):
      SlaterDetOperations_base<AllocType>(NMO,NAEA,alloc_)
    {
    }

    ~SlaterDetOperations_serial() {}

    SlaterDetOperations_serial(const SlaterDetOperations_serial& other) = delete;
    SlaterDetOperations_serial(SlaterDetOperations_serial&& other) = default;
    SlaterDetOperations_serial& operator=(const SlaterDetOperations_serial& other) = delete;
    SlaterDetOperations_serial& operator=(SlaterDetOperations_serial&& other) = default;

    // C must live in shared memory for this routine to work as expected
    template<class MatA, class MatB, class MatC>
    void MixedDensityMatrix(const MatA& hermA, const MatB& B, MatC&& C, T* res, communicator& comm, bool compact=false) {
#ifdef QMC_CUDA
      APP_ABORT(" Error: SlaterDetOperations_serial should not be here. \n");
#endif
      Base::MixedDensityMatrix(hermA,B,C,res,compact);
    }

    template<class integer, class MatA, class MatB, class MatC, class MatQ>
    void MixedDensityMatrixForWoodbury(const MatA& hermA, const MatB& B, MatC&& C, T* res,
                                    integer* ref, MatQ&& QQ0, communicator& comm, bool compact=false) {
#ifdef QMC_CUDA
      APP_ABORT(" Error: SlaterDetOperations_serial should not be here. \n");
#endif
      Base::MixedDensityMatrixForWoodbury(hermA,B,std::forward<MatC>(C),res,ref,std::forward<MatQ>(QQ0),
                                    compact);
    }

    template<class MatA, class MatB>
    void Overlap(const MatA& hermA, const MatB& B, T* res, communicator& comm) { 
#ifdef QMC_CUDA
      APP_ABORT(" Error: SlaterDetOperations_serial should not be here. \n");
#endif
      Base::Overlap(hermA,B,res);
    }

    template<typename integer, class MatA, class MatB, class MatC>
    void OverlapForWoodbury(const MatA& hermA, const MatB& B, T* res, integer* ref, MatC&& QQ0, communicator& comm) {
#ifdef QMC_CUDA
      APP_ABORT(" Error: SlaterDetOperations_serial should not be here. \n");
#endif
      Base::OverlapForWoodbury(hermA,B,res,ref,std::forward<MatC>(QQ0));
    }

    template<class Mat, class MatP1, class MatV>
    void Propagate(Mat&& A, const MatP1& P1, const MatV& V, communicator& comm, int order=6) {
#ifdef QMC_CUDA
      APP_ABORT(" Error: SlaterDetOperations_serial should not be here. \n");
#endif
      Base::Propagate(std::forward<Mat>(A),P1,V,order);
    }

    template<class MatA, class MatP1, class MatV>
    void BatchedPropagate(std::vector<MatA> &Ai, const MatP1& P1, const MatV& V, int order=6) {
      static_assert( std::decay<MatA>::type::dimensionality == 2, " dimenionality == 2" );
      static_assert( std::decay<MatV>::type::dimensionality == 3, " dimenionality == 3" );
      if(Ai.size() == 0) return;
      assert(Ai.size() == V.size(0));
      int nbatch = Ai.size();
      int NMO = Ai[0].size(0);
      int NAEA = Ai[0].size(1);
//      for(int i=0; i<Ai.size(); ++i)
//        Propagate(Ai[i],P1,V[i],order);
      set_buffer(nbatch*3*NMO*NAEA);
      TTensor_ref TMN(TBuff.data(), {nbatch,NMO,NAEA});
      TTensor_ref T1(TMN.data()+TMN.num_elements(), {nbatch,NMO,NAEA});
      TTensor_ref T2(T1.data()+T1.num_elements(), {nbatch,NMO,NAEA});
      // could be batched when csrmm is batched  
      for(int ib=0; ib<nbatch; ib++) 
        ma::product(P1,Ai[ib],TMN[ib]);
      SlaterDeterminantOperations::batched::apply_expM(V,TMN,T1,T2,order);
      for(int ib=0; ib<nbatch; ib++) 
        ma::product(P1,TMN[ib],Ai[ib]);
    }

    // C[nwalk, M, N]
    template<class MatA, class MatB, class MatC, class TVec>
    void BatchedMixedDensityMatrix(const MatA& hermA, std::vector<MatB> &Bi, MatC&& C, TVec&& ovlp, bool compact=false) {
      if(Bi.size()==0) return;
      static_assert(std::decay<MatC>::type::dimensionality == 3, "Wrong dimensionality");
      static_assert(std::decay<TVec>::type::dimensionality == 1, "Wrong dimensionality");
      int NMO = hermA.size(1);
      int NAEA = hermA.size(0);
      int nbatch = Bi.size();
      assert(C.size(0) == nbatch);
      assert(ovlp.size() == nbatch);
      int n1=nbatch, n2=NAEA, n3=NMO;
      if(compact) {
        n1=n2=n3=0;
      }
      set_buffer(nbatch*NAEA*NAEA + n1*n2*n3);
      TTensor_ref TNN3D(TBuff.origin(), {nbatch,NAEA,NAEA});
      TTensor_ref TNM3D(TBuff.origin()+TNN3D.num_elements(), {n1,n2,n3});
      int sz = ma::invert_optimal_workspace_size(TNN3D[0]);
      if(WORK.num_elements() < nbatch*sz)
        WORK.reextent(iextensions<1u>{nbatch*sz});
      if(IWORK.num_elements() < nbatch*(NMO+1))
        IWORK.reextent(iextensions<1u>{nbatch*(NMO+1)});
      SlaterDeterminantOperations::batched::MixedDensityMatrix(hermA,Bi,
                std::forward<MatC>(C),std::forward<TVec>(ovlp),TNN3D,TNM3D,IWORK,WORK,compact);
    }

    template<class MatA, class MatB, class TVec>
    void BatchedOverlap(const MatA& hermA, std::vector<MatB> &Bi, TVec&& ovlp) {
      if(Bi.size()==0) return;
      static_assert(std::decay<TVec>::type::dimensionality == 1, "Wrong dimensionality");
      int NMO = hermA.size(1);
      int NAEA = hermA.size(0);
      int nbatch = Bi.size();
      assert(ovlp.size() == nbatch);
      set_buffer(nbatch*NAEA*NAEA);
      TTensor_ref TNN3D(TBuff.origin(), {nbatch,NAEA,NAEA});
      if(IWORK.num_elements() < nbatch*(NMO+1))
        IWORK.reextent(iextensions<1u>{nbatch*(NMO+1)});
      SlaterDeterminantOperations::batched::Overlap(hermA,Bi,
                std::forward<TVec>(ovlp),TNN3D,IWORK);
    }

    template<class MatA>
    void BatchedOrthogonalize(std::vector<MatA> &Ai, T* detR=nullptr) {
/*
      if(detR == nullptr)
        for(int i=0; i<Ai.size(); i++)
          Orthogonalize(Ai[i],detR);
      else  
        for(int i=0; i<Ai.size(); i++)
          Orthogonalize(Ai[i],detR+i);
      return;
*/
#ifdef QMC_CUDA
      // QR on the transpose
      if(Ai.size()==0) return;
      int NMO = Ai[0].size(0);
      int NAEA = Ai[0].size(1);
      int nbatch = Ai.size();
      set_buffer(nbatch*(NMO*NAEA + 2*NMO));
      TTensor_ref AT(TBuff.origin(), {nbatch,NAEA,NMO});
      TMatrix_ref T_(AT.origin() + AT.num_elements(), {nbatch,NMO});
      TMatrix_ref scl(T_.origin() + T_.num_elements(), {nbatch,NMO});
      if(IWORK.num_elements() < nbatch*(NMO+1))
        IWORK.reextent(iextensions<1u>{nbatch*(NMO+1)});
      int sz = ma::gqr_optimal_workspace_size(AT[0]);
      if(WORK.num_elements() < nbatch*sz)
        WORK.reextent(iextensions<1u>{nbatch*sz});
      for(int i=0; i<nbatch; i++)
        ma::transpose(Ai[i],AT[i]);
      // careful, expects fortran order
      geqrfStrided(NMO,NAEA,AT.origin(),NMO,NMO*NAEA,T_.origin(),NMO,IWORK.origin(),nbatch);
      //for(int i=0; i<nbatch; i++)
      //  ma::geqrf(AT[i],T_[i],WORK); 
      using ma::determinant_from_geqrf;
      using ma::scale_columns;
// needs batched
      for(int i=0; i<nbatch; i++) {
        if(detR != nullptr)
          determinant_from_geqrf(NAEA,AT[i].origin(),NMO,scl[i].origin(),detR+i);
        else 
          determinant_from_geqrf(NAEA,AT[i].origin(),NMO,scl[i].origin());
      }
//      for(int i=0; i<nbatch; i++) 
//        ma::gqr(AT[i],T_[i],WORK);
      gqrStrided(NMO,NAEA,NAEA,AT.origin(),NMO,NMO*NAEA,T_.origin(),NMO,WORK.origin(),sz,IWORK.origin(),nbatch);
      for(int i=0; i<nbatch; i++) {
        ma::transpose(AT[i],Ai[i]);
        scale_columns(NMO,NAEA,Ai[i].origin(),Ai[i].stride(0),scl[i].origin());
      }
#else
      int nw=Ai.size();
      for(int i=0; i<nw; i++) 
        Orthogonalize(Ai[i], detR+i);
#endif
    }

  protected:

    using Base::IWORK;
    using Base::WORK;
    using Base::TBuff;
    using Base::set_buffer;

};

}

}

#endif
