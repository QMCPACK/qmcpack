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

#ifndef QMCPLUSPLUS_AFQMC_SERIALPROPAGATOR_H
#define QMCPLUSPLUS_AFQMC_SERIALPROPAGATOR_H

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

#include "AFQMC/Wavefunctions/Wavefunction_serial.hpp"

namespace qmcplusplus
{

namespace afqmc
{

/*
 * AFQMC propagator using shared memory. Data structures are in shared memory, 
 * but not distributed over multiple nodes. 
 */
template<class Alloc>
class AFQMCSerialPropagator: public AFQMCInfo
{
  protected:

  using Allocator = Alloc; 
  using T = typename Allocator::value_type;
  using pointer = typename Allocator::pointer;
  using const_pointer = typename Allocator::const_pointer;

  using CVector = boost::multi::array<ComplexType,1,Allocator>;  
  using CVector_ref = boost::multi::array_ref<ComplexType,1,pointer>;  
  using CMatrix = boost::multi::array<ComplexType,2,Allocator>;  
  using CMatrix_ref = boost::multi::array_ref<ComplexType,2,pointer>;  
  using C3Tensor = boost::multi::array<ComplexType,3,Allocator>;  
  using C3Tensor_ref = boost::multi::array_ref<ComplexType,3,pointer>;  

  public:

    AFQMCSerialPropagator(AFQMCInfo& info, xmlNodePtr cur, afqmc::TaskGroup_& tg_, 
                          Wavefunction_serial& wfn_, CMatrix&& h1_, CVector&& vmf_, 
                          RandomGenerator_t* r, Allocator alloc_={}):
            AFQMCInfo(info),allocator_(alloc_),TG(tg_),wfn(wfn_),
            H1(std::move(h1_)),
            //P1(P1Type(tp_ul_ul{0,0},tp_ul_ul{0,0},0,boost::mpi3::intranode::allocator<ComplexType>(tg_.TG_local()))),
            P1({NMO,NMO},allocator_),
            rng(r),
            //SDetOp(2*NMO,NAEA+NAEB),
            SDetOp(SlaterDetOperations_serial<Allocator>(2*NMO,NAEA+NAEB,allocator_)),
            buff(extensions<1u>{1},allocator_),
            local_vHS({0,0},allocator_),
            new_overlaps(extensions<1u>{0},allocator_), 
            new_energies({0,0},allocator_), 
            MFfactor({0,0},allocator_), 
            hybrid_weight({0,0},allocator_), 
            vMF(std::move(vmf_)),
            TSM({2*NMO,NAEA+NAEB},allocator_), 
            old_dt(-123456.789),
            order(6),
            nbatched_propagation(-1)
    {
      transposed_vHS_ = wfn.transposed_vHS();
      transposed_G_ = wfn.transposed_G_for_vbias();
      if(not transposed_vHS_) local_vHS.reextent({NMO,NMO});
      parse(cur);  
    }

    ~AFQMCSerialPropagator() {}

    AFQMCSerialPropagator(AFQMCSerialPropagator const& other) = delete;
    AFQMCSerialPropagator& operator=(AFQMCSerialPropagator const& other) = delete;
    AFQMCSerialPropagator(AFQMCSerialPropagator&& other) = default;
    AFQMCSerialPropagator& operator=(AFQMCSerialPropagator&& other) = default;

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

    // reset shared memory buffers
    // useful when the current buffers use too much memory (e.g. reducing steps in future calls)
    void reset() { shmbuff.reextent(extensions<1u>{0}); }

    bool hybrid_propagation() { return hybrid; }

    bool free_propagation() { return free_projection; }

  protected: 

    Allocator allocator_;

    TaskGroup_& TG;

    Wavefunction_serial& wfn;

    // P1 = exp(-0.5*dt*H1), so H1 includes terms from MF substraction 
    //                       and the exchange term from the cholesky decomposition (e.g. vn0) 
    CMatrix H1;

    P1Type P1;

    RandomGenerator_t* rng;

    SlaterDetOperations_shared<ComplexType> SDetOp;
    //SlaterDetOperations SDetOp;

    CVector buff;    

    RealType old_dt;
    int order;

    RealType vbias_bound;

    // type of propagation
    int nbatched_propagation;
    bool free_projection;
    bool hybrid;
    bool importance_sampling;
    bool apply_constrain;

    bool transposed_vHS_;
    bool transposed_G_;

    // used in propagator step 
    CMatrix local_vHS;  

    CVector new_overlaps;  
    CMatrix new_energies;  

    CMatrix MFfactor;  
    CMatrix hybrid_weight;  

    CVector vMF;  
    // Temporary for propagating with constructed B matrix.
    CMatrix TSM;
 
    template<class WlkSet>
    void step(int steps, WlkSet& wset, RealType E1, RealType dt);

    void assemble_X(size_t nsteps, size_t nwalk, RealType sqrtdt,
                    CMatrix_ref& X, CMatrix_ref & vbias,
                    CMatrix_ref& MF, CMatrix_ref& HW, bool addRAND=true);

    void parse(xmlNodePtr cur);

    template<class WSet>
    void apply_propagators(WSet& wset, int ni, C3Tensor_ref& vHS3D); 

    template<class WSet>
    void apply_propagators_construct_propagator(WSet& wset, int ni, C3Tensor_ref& vHS3D); 

    ComplexType apply_bound_vbias(ComplexType v, RealType sqrtdt)
    {
// Needs kernel!!!
      return (std::abs(v)>vbias_bound*sqrtdt)?
                (v/(std::abs(v)/static_cast<ValueType>(vbias_bound*sqrtdt))):(v);
    }
  
};

}

}

#include "AFQMCSerialPropagator.icc"

#endif

