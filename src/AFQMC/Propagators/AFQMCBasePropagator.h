//////////////////////////////////////////////////////////////////////
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

#ifndef QMCPLUSPLUS_AFQMC_BASEPROPAGATOR_H
#define QMCPLUSPLUS_AFQMC_BASEPROPAGATOR_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <tuple>

#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/RandomGenerator.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/Propagators/generate1BodyPropagator.hpp"

#include "AFQMC/Wavefunctions/Wavefunction.hpp"

namespace qmcplusplus
{

namespace afqmc
{

/*
 * Base class for AFQMC propagators.
 */
class AFQMCBasePropagator: public AFQMCInfo
{
  protected:

  // allocator for local memory
  using allocator = device_allocator<ComplexType>;
  using pointer = device_ptr<ComplexType>;
  // allocator for memory shared by all cores in working local group 
  using aux_allocator = localTG_allocator<ComplexType>;

  using CVector = boost::multi::array<ComplexType,1,allocator>;  
  using CMatrix = boost::multi::array<ComplexType,2,allocator>;  
  using C3Tensor = boost::multi::array<ComplexType,3,allocator>;  
  using CVector_ref = boost::multi::array_ref<ComplexType,1,pointer>;  
  using CMatrix_ref = boost::multi::array_ref<ComplexType,2,pointer>;  
  using C3Tensor_ref = boost::multi::array_ref<ComplexType,3,pointer>;  
  using sharedCVector = ComplexVector<aux_allocator>; 
  using stdCVector = boost::multi::array<ComplexType,1>;  
  using stdCMatrix = boost::multi::array<ComplexType,2>;  
  using stdCVector_ref = boost::multi::array_ref<ComplexType,1>;  
  using stdCMatrix_ref = boost::multi::array_ref<ComplexType,2>;  
  using stdC3Tensor_ref = boost::multi::array_ref<ComplexType,3>;  

  using mpi3CVector = boost::multi::array<ComplexType,1,shared_allocator<ComplexType>>;  
  using mpi3CMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;  
  using mpi3CTensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;  

  public:

    AFQMCBasePropagator(AFQMCInfo& info, xmlNodePtr cur, afqmc::TaskGroup_& tg_, 
                          Wavefunction& wfn_, 
                          stdCMatrix&& h1_, CVector&& vmf_, 
                          RandomGenerator_t* r):
            AFQMCInfo(info),TG(tg_),
            alloc_(),aux_alloc_(make_localTG_allocator<ComplexType>(TG)),wfn(wfn_),
            H1(std::move(h1_),shared_allocator<ComplexType>{TG.Node()}),
            H1ext({2,1,1},shared_allocator<ComplexType>{TG.Node()}),
            vMF(std::move(vmf_)),
            rng(r),
            SDetOp(wfn.getSlaterDetOperations()),
            buffer(iextensions<1u>{1},aux_alloc_),
            local_group_comm(),
            old_dt(-123456.789),
            last_nextra(-1),
            last_task_index(-1),
            order(6),
            nbatched_propagation(0),
            nbatched_qr(0),
            spin_dependent_P1(false)
    {
      P1.reserve(2);  
      P1.emplace_back(P1Type(tp_ul_ul{0,0},tp_ul_ul{0,0},0,aux_alloc_));
      transposed_vHS_ = wfn.transposed_vHS();
      transposed_G_ = wfn.transposed_G_for_vbias();
      parse(cur);  
    }

    ~AFQMCBasePropagator() {}

    AFQMCBasePropagator(AFQMCBasePropagator const& other) = delete;
    AFQMCBasePropagator& operator=(AFQMCBasePropagator const& other) = delete;
    AFQMCBasePropagator(AFQMCBasePropagator&& other) = default;
    AFQMCBasePropagator& operator=(AFQMCBasePropagator&& other) = default;

    template<class WlkSet>
    void Propagate(int steps, WlkSet& wset, RealType E1,
                   RealType dt, int fix_bias=1) { 
      int nblk = steps/fix_bias;
      int nextra = steps%fix_bias;
      for(int i=0; i<nblk; i++) {
        step(fix_bias,wset,E1,dt);
      }
      if(nextra>0) 
        step(nextra,wset,E1,dt);
      TG.local_barrier();
    }

    template<class WlkSet, class CTens, class CMat>
    void BackPropagate(int steps, int nStabalize, WlkSet& wset, CTens&& Refs, CMat&& detR); 

    // reset shared memory buffers
    // useful when the current buffers use too much memory (e.g. reducing steps in future calls)
    void reset() { buffer.reextent(iextensions<1u>{0}); }

    bool hybrid_propagation() { return hybrid; }

    bool free_propagation() { return free_projection; }

    int global_number_of_cholesky_vectors() const{
      return wfn.global_number_of_cholesky_vectors(); 
    }

    // in case P1 needs to exist before call to Propagate is executed
    void generateP1(double dt, WALKER_TYPES walker_type) {
      old_dt = dt;
      // generate1BodyPropagator currently expects a shared_allocator, fix later
      using P1shm = ma::sparse::csr_matrix<ComplexType,int,int,
                                shared_allocator<ComplexType>,
                                ma::sparse::is_root>;
      if(spin_dependent_P1) {
        if(walker_type!=COLLINEAR)
           APP_ABORT(" Error: Spin dependent P1 being used with CLOSED walker.\n");
        P1[0] = std::move(generate1BodyPropagator<P1shm>(TG,1e-8,dt,H1,H1ext[0],printP1eV));
        P1[1] = std::move(generate1BodyPropagator<P1shm>(TG,1e-8,dt,H1,H1ext[1],printP1eV));
      } else {
        P1[0] = std::move(generate1BodyPropagator<P1shm>(TG,1e-8,dt,H1,printP1eV));
      } 
    }

  protected: 

    TaskGroup_& TG;

    allocator alloc_;
 
    aux_allocator aux_alloc_;

    Wavefunction& wfn;

    // P1 = exp(-0.5*dt*H1), so H1 includes terms from MF substraction 
    //                       and the exchange term from the cholesky decomposition (e.g. vn0) 
    mpi3CMatrix H1;
    mpi3CTensor H1ext;

    std::vector<P1Type> P1;

    CVector vMF;  

    RandomGenerator_t* rng;

    SlaterDetOperations* SDetOp;

    sharedCVector buffer;    

    shared_communicator local_group_comm;

    RealType old_dt;
    int last_nextra;
    int last_task_index;
    int order;
    int nbatched_propagation;
    int nbatched_qr;
    bool spin_dependent_P1;
    bool printP1eV=false;

    RealType vbias_bound;

    // type of propagation
    bool free_projection;
    bool hybrid;
    bool importance_sampling;
    bool apply_constrain;

    bool transposed_vHS_;
    bool transposed_G_;

    // used in propagator step 
    // intead of doing this, should use TBuff to transpose vHS3D and only have propagation 
    // with vHD3D[nwalk*nsteps,...]
    CMatrix local_vHS;  

    CVector new_overlaps;  
    CMatrix new_energies;  

    CMatrix MFfactor;  
    CMatrix hybrid_weight;  

    boost::multi::array<ComplexType,2> work; 
 
    template<class WlkSet>
    void step(int steps, WlkSet& wset, RealType E1, RealType dt);

    void assemble_X(size_t nsteps, size_t nwalk, RealType sqrtdt,
                    CMatrix_ref& X, CMatrix_ref & vbias,
                    CMatrix_ref& MF, CMatrix_ref& HW, bool addRAND=true);

    void reset_nextra(int nextra); 

    void parse(xmlNodePtr cur);

    template<class WSet>
    void apply_propagators(char TA, WSet& wset, int ni, int tk0, int tkN, int ntask_total_serial,
                           C3Tensor_ref& vHS3D);

    template<class WSet>
    void apply_propagators_batched(char TA, WSet& wset, int ni, C3Tensor_ref& vHS3D);

    ComplexType apply_bound_vbias(ComplexType v, RealType sqrtdt)
    {
      return (std::abs(v)>vbias_bound*sqrtdt)?
                (v/(std::abs(v)/static_cast<ValueType>(vbias_bound*sqrtdt))):(v);
    }

    // taken from NOMSD 
    template<class WlkSet, class CMat>
    void Orthogonalize_batched(WlkSet& wset, CMat&& detR);

    template<class WlkSet, class CMat>
    void Orthogonalize_shared(WlkSet& wset, CMat&& detR);

};

}

}

#include "AFQMCBasePropagator.icc"

#endif

