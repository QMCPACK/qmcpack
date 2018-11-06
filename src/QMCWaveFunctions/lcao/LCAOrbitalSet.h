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
    ValueMatrix_t* m_init_B;
    ValueMatrix_t  T;
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
    ///test
    void test()
    {
      app_log()<<"test function test function SDP \n";
    }
    /// create optimizable orbital rotation parameters
    void buildOptVariables(std::vector<RealType>& input_params, bool params_supplied, std::vector<int> * data, const size_t& nel, std::vector<size_t>& C2node, const int& spin);

    ///helper function to buildOptVariables
    int build_occ_vec(std::vector<int> * data, const size_t& nel, const size_t& nmo, std::vector<int>* occ_vec);

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
      app_log() << "RESET PARAMETERS CALLED SDP\n";
      //myBasisSet->resetParameters(active);
      if (Optimizable)
      {
        const size_t& nmo = OrbitalSetSize;
        const size_t& nb = BasisSetSize;
        // read out the parameters for spin up electrons that define the rotation into an antisymmetric matrix
        std::vector<RealType> rot_mat(nmo*nmo, 0.0);
        for (int i=0; i < m_act_rot_inds.size(); i++)
        {
          int loc=myVars.where(i);

          const int p = m_act_rot_inds[i].first;
          const int q = m_act_rot_inds[i].second;
          // m_first_var_pos is the index to the first parameter of the spin up electrons...
          const RealType x = myVars[i] = active[loc];
          //app_log() <<"active ["<< loc << "] = "<< x<< "\n";
          rot_mat[p+q*nmo] =  x;
          rot_mat[q+p*nmo] = -x;
        }

        this->exponentiate_antisym_matrix(nmo, &rot_mat[0]);

        BLAS::gemm('N','T', nb, nmo, nmo, RealType(1.0), m_init_B->data(), nb, &rot_mat[0], nmo, RealType(0.0), C->data(), nb);

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
    //helper function to resetParameters
    void exponentiate_antisym_matrix(const int n, RealType* const mat);

  };
}
#endif
