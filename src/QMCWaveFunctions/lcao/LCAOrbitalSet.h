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

#include <Numerics/MatrixOperators.h>
#include "Numerics/DeterminantOperators.h"

namespace qmcplusplus
{
  /** class to handle linear combinations of basis orbitals used to evaluate the Dirac determinants.
   *
   * SoA verson of LCOrtbitalSet
   * Localized basis set is always real 
   */
  struct LCAOrbitalSet: public SPOSet
  {
  public:
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
    // The initial coefficents at the start of the simulation
    ValueMatrix_t* m_init_B;

    ///true if C is an identity matrix
    bool Identity;
    ///if true, do not clean up
    bool IsCloned;
    ///Temp(BasisSetSize) : Row index=V,Gx,Gy,Gz,L
    vgl_type Temp; 
    ///Tempv(OrbitalSetSize) Tempv=C*Temp
    vgl_type Tempv; 
    //lookup table mapping the unique determinants to their element position in C2_node vector
    std::vector< std::vector<int> > lookup_tbl;
    //vector that contains active orbital rotation parameter indices 
    std::vector<std::pair<int,int> > m_act_rot_inds;
    /** constructor
     * @param bs pointer to the BasisSet
     * @param rl report level
     */
    LCAOrbitalSet(basis_type* bs=nullptr,int rl=0);

    LCAOrbitalSet(const LCAOrbitalSet& in)=default;

    virtual ~LCAOrbitalSet();

    SPOSet* makeClone() const;

    /// create optimizable orbital rotation parameters
    void buildOptVariables(std::vector<int> * data, const size_t& nel, std::vector<size_t>& C2node, const int& spin);

    ///helper function to buildOptVariables
    int build_occ_vec(std::vector<int> * data, const size_t& nel, const size_t& nmo, std::vector<int>* occ_vec);
  
    void evaluateDerivatives (ParticleSet& P,
                             const opt_variables_type& optvars,
                             std::vector<RealType>& dlogpsi, 
                             std::vector<RealType>& dhpsioverpsi,
                             const ValueType& psiCurrent,
                             std::vector<RealType> const * const Coeff,
                             std::vector<size_t> const * const C2node_up,
                             std::vector<size_t> const * const C2node_dn,
                             const ValueVector_t& detValues_up, 
                             const ValueVector_t& detValues_dn, 
                             const GradMatrix_t& grads_up, 
                             const GradMatrix_t& grads_dn, 
                             const ValueMatrix_t& lapls_up, 
                             const ValueMatrix_t& lapls_dn,
                             const ValueMatrix_t& M_up,
                             const ValueMatrix_t& M_dn,
                             const ValueMatrix_t& Minv_up,
                             const ValueMatrix_t& Minv_dn,
                             const GradMatrix_t& B_grad,
                             const ValueMatrix_t& B_lapl,
                             std::vector<int> const * const detData_up, 
                             const size_t N1,
                             const size_t N2,
                             const size_t NP1,
                             const size_t NP2);

    
    void checkInVariables(opt_variables_type& active)
    {
      if(Optimizable && !IsCloned)
      {
        if(myVars.size())
          active.insertFrom(myVars);
        else
          Optimizable=false;
      }
    }

    void checkOutVariables(const opt_variables_type& active)
    {
      if(Optimizable && !IsCloned)
        myVars.getIndex(active);
    }


    ///reset
    void resetParameters(const opt_variables_type& active)
    {
      if (Optimizable)
      {
        std::vector<RealType> param( m_act_rot_inds.size() );
        for (int i=0; i < m_act_rot_inds.size(); i++)
        {
          int loc=myVars.where(i);
          param[i] = myVars[i] = active[loc];
        }
        apply_rotation(&param);

      }

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

  private:
    //fucntion to perform orbital rotations
    void apply_rotation(std::vector<RealType> const * const param);

    //helper function to apply_rotation 
    void exponentiate_antisym_matrix(const int n, RealType* const mat);
    

    //helper function to evaluatederivative; evaluate orbital rotation parameter derivative using table method
    void table_method_eval(std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi,
                           const ParticleSet::ParticleLaplacian_t& myL_J,
                           const ParticleSet::ParticleGradient_t& myG_J,
                           const size_t nel,
                           const size_t nmo,
                           const ValueType& psiCurrent,
                           std::vector<RealType> const * const Coeff,
                           std::vector<size_t> const * const C2node_up,
                           std::vector<size_t> const * const C2node_dn,
                           const ValueVector_t& detValues_up, 
                           const ValueVector_t& detValues_dn, 
                           const GradMatrix_t& grads_up, 
                           const GradMatrix_t& grads_dn, 
                           const ValueMatrix_t& lapls_up, 
                           const ValueMatrix_t& lapls_dn,
                           const ValueMatrix_t& M_up,
                           const ValueMatrix_t& M_dn,
                           const ValueMatrix_t& Minv_up,
                           const ValueMatrix_t& Minv_dn,
                           const GradMatrix_t& B_grad,
                           const ValueMatrix_t& B_lapl,
                           std::vector<int> const * const detData_up,
                           const size_t N1,
                           const size_t N2,
                           const size_t NP1,
                           const size_t NP2); 

  };
}
#endif
