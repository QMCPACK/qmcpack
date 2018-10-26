//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_SOA_LINEARCOMIBINATIONORBITALSET_TEMP_H
#define QMCPLUSPLUS_SOA_LINEARCOMIBINATIONORBITALSET_TEMP_H

#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/BasisSetBase.h"


namespace qmcplusplus
{
  /** class to handle linear combinations of basis orbitals used to evaluate the Dirac determinants.
   *
   * SoA verson of LCOrtbitalSet
   * Localized basis set is always real 
   */
  struct LCAOrbitalSet: public SPOSet
  {
    typedef RealBasisSetBase<RealType> basis_type;
    typedef basis_type::vgl_type vgl_type;

    ///level of printing
    int ReportLevel;
    ///pointer to the basis set
    basis_type* myBasisSet;
    ///number of Single-particle orbitals
    IndexType BasisSetSize;
    /** pointer to matrix containing the coefficients
     *
     * makeClone makes a shallow copy
     */
    ValueMatrix_t* C;
    ///true if C is an identity matrix
    bool Identity;
    ///if true, do not clean up
    bool IsCloned;
    ///Temp(BasisSetSize) : Row index=V,Gx,Gy,Gz,L
    vgl_type Temp; 
    ///Tempv(OrbitalSetSize) Tempv=C*Temp
    vgl_type Tempv; 

    /** constructor
     * @param bs pointer to the BasisSet
     * @param rl report level
     */
    LCAOrbitalSet(basis_type* bs=nullptr,int rl=0);

    LCAOrbitalSet(const LCAOrbitalSet& in)=default;

    virtual ~LCAOrbitalSet();

    SPOSet* makeClone() const;

    ///reset
    void resetParameters(const opt_variables_type& active)
    {
      //myBasisSet->resetParameters(active);
    }

    ///reset the target particleset
    void resetTargetParticleSet(ParticleSet& P)
    {
      //myBasisSet->resetTargetParticleSet(P);
    }

    /** set the OrbitalSetSize
    */
    virtual void setOrbitalSetSize(int norbs)
    {
      OrbitalSetSize=norbs;
      Tempv.resize(OrbitalSetSize);
    }

    /** set the basis set
    */
    void setBasisSet(basis_type* bs);

    /** return the size of the basis set
    */
    int getBasisSetSize() const
    {
      return (myBasisSet==nullptr)? 0: myBasisSet->getBasisSetSize();
    }

    bool setIdentity(bool useIdentity);

    void checkObject() const
    {
      if(!(OrbitalSetSize == C->rows() && BasisSetSize == C->cols()))
        APP_ABORT("   LCAOrbitalSet::checkObject Linear coeffient for LCAOrbitalSet is not consistent with the input.");
    }

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

    void evaluateValues(const VirtualParticleSet& VP, ValueMatrix_t& psiM, ValueAlignedVector_t& SPOMem);

    size_t estimateMemory(const int nP);

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi);

    void evaluate_notranspose(const ParticleSet& P, int first, int last,
        ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);

    void evaluate_notranspose(const ParticleSet& P, int first, int last,
        ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet);

    void evaluate_notranspose(const ParticleSet& P, int first, int last
        , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet);

    void evaluateThirdDeriv(const ParticleSet& P, int first, int last , GGGMatrix_t& grad_grad_grad_logdet);

    //helper functions to handl Identity
    void evaluate_vgl_impl(const vgl_type& temp,
        ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi) const;
    void evaluate_vgl_impl(const vgl_type& temp, int i,
        ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet) const;

  };
}
#endif
