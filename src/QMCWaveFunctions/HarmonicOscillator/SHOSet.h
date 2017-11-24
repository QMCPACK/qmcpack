//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef QMCPLUSPLUS_SHOSET_H
#define QMCPLUSPLUS_SHOSET_H

#include <QMCWaveFunctions/SPOSetBase.h>
#include <QMCWaveFunctions/SPOInfo.h>
#include <Numerics/VectorViewer.h>


namespace qmcplusplus
{

  struct SHOState : public SPOInfo
  {

    TinyVector<int,DIM> quantum_number;

    SHOState()
    {
      quantum_number = -1;
      energy         = 0.0;
    }

    ~SHOState() { }

    inline void set(TinyVector<int,DIM> qn, RealType e)
    {
      quantum_number = qn;
      energy         = e;
    }

    inline void sho_report(const std::string& pad="")
    {
      app_log()<<pad<<"qn="<<quantum_number<<"  e="<<energy<< std::endl;
    }
  };


  struct SHOSet : public SPOSetBase
  {
    typedef ValueMatrix_t::value_type value_type;
    typedef GradMatrix_t::value_type grad_type;
    typedef VectorViewer<value_type> ValueViewer_t;
    typedef VectorViewer<grad_type>  GradViewer_t;

    RealType length; 
    PosType  center;

    int nmax;
    TinyVector<int,DIM> qn_max;
    std::vector<SHOState>  state_info;
    std::vector<RealType>  prefactors;
    Array<RealType,2> hermite;
    Array<RealType,2> bvalues;
    Array<RealType,2> d0_values;
    Array<RealType,2> d1_values;
    Array<RealType,2> d2_values;

    //construction/destruction
    SHOSet(RealType l, PosType c, const std::vector<SHOState*>& sho_states);

    ~SHOSet();

    void initialize();


    //SPOSetBase interface methods
    SPOSetBase* makeClone() const;

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi,
                  GradVector_t& dpsi, ValueVector_t& d2psi);

    void evaluate_notranspose(const ParticleSet& P, int first, int last, 
                              ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
                              ValueMatrix_t& d2logdet);


    //local functions
    void evaluate_v(PosType r, ValueViewer_t& psi);
    void evaluate_vgl(PosType r, ValueViewer_t& psi, 
                      GradViewer_t& dpsi, ValueViewer_t& d2psi);
    void evaluate_hermite(const PosType& xpos);
    void evaluate_d0(const PosType& xpos, ValueViewer_t& psi);
    void evaluate_d1(const PosType& xpos, ValueViewer_t& psi, GradViewer_t& dpsi);
    void evaluate_d2(const PosType& xpos, ValueViewer_t& psi, ValueViewer_t& d2psi);
    void report(const std::string& pad="");
    void test_derivatives();
    void test_overlap();
    void evaluate_check(PosType r, ValueViewer_t& psi, 
                        GradViewer_t& dpsi, ValueViewer_t& d2psi);

    //empty methods 
    /// number of orbitals is determined only by initial request
    inline void setOrbitalSetSize(int norbs)                     { }

    /// does not affect ParticleSet information
    inline void resetTargetParticleSet(ParticleSet& P)           { }

    
    ///unimplemented functions call this to abort
    inline void not_implemented(const std::string& method)
    {
      APP_ABORT("SHOSet::"+method+" has not been implemented.");
    }


    //methods to be implemented in the future (possibly)
    void resetParameters(const opt_variables_type& optVariables);
    void evaluate(const ParticleSet& P, PosType &r, ValueVector_t &psi);
    void evaluate(const ParticleSet& P, int iat,ValueVector_t& psi, 
                  GradVector_t& dpsi, HessVector_t& ddpsi);
    void evaluateThirdDeriv(const ParticleSet& P,int first,int last,
                            GGGMatrix_t& dddlogdet);
    void evaluate_notranspose(const ParticleSet& P, int first, int last, 
                              ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
                              HessMatrix_t& ddlogdet);
    void evaluate_notranspose(const ParticleSet& P, int first, int last, 
                              ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
                              HessMatrix_t& ddlogdet, GGGMatrix_t& dddlogdet);
    void evaluateGradSource(const ParticleSet &P, int first, int last, 
                            const ParticleSet &source, int iat_src, 
                            GradMatrix_t &gradphi);
    void evaluateGradSource(const ParticleSet &P, int first, int last, 
                            const ParticleSet &source, int iat_src, 
                            GradMatrix_t &dphi, HessMatrix_t &ddphi, 
                            GradMatrix_t &dlapl_phi);
    void evaluateBasis(const ParticleSet &P, int first, int last, 
                       ValueMatrix_t &basis_val, GradMatrix_t &basis_grad, 
                       ValueMatrix_t &basis_lapl);
    void evaluateForDeriv(const ParticleSet &P, int first, int last, 
                          ValueMatrix_t &basis_val, GradMatrix_t &basis_grad, 
                          ValueMatrix_t &basis_lapl);
    void copyParamsFromMatrix(const opt_variables_type& active, 
                              const ValueMatrix_t &mat, 
                              std::vector<RealType> &destVec);
  };

}


#endif
