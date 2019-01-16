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
    void buildOptVariables(std::vector<RealType>& input_params, bool params_supplied, std::vector<int> * data, const size_t& nel, std::vector<size_t>& C2node, const int& spin);

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

        exponentiate_antisym_matrix(nmo, &rot_mat[0]);

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

  private:
    //helper function to resetParameters
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

    //helper function to contruct tables 'T' and 'M', which must be computed for both the single slater and multi-slater wfns.
    //Also computes the special O operator matrix that is needed in both cases as well. 
    void construct_tables(
                          const ParticleSet::ParticleLaplacian_t& myL_J,
                          const ParticleSet::ParticleGradient_t& myG_J,
                          const GradMatrix_t& B_grad,
                          const ValueMatrix_t& B_lapl,
                          const ValueMatrix_t& M_up,
                          const ValueMatrix_t& Minv_up,
                          const size_t nel,
                          const size_t nmo,
                          const int offset1,
                          double* T
                          );

    /*
    ##############MATRICES DEFINED
    ##############The following is also presented in the table_method_eval function in the context
    ##############of explaining how the table method works.
    Below is a translation of the shorthand I use to represent matrices independent of ``excitation matrix".
    \newline
    \hfill\break
    $
        Y_{1} =  A^{-1}B   \\
        Y_{2} = A^{-1}BA^{-1}\widetilde{A} \\
        Y_{3} = A^{-1}\widetilde{B} \\
        Y_{4} = \widetilde{M} = (A^{-1}\widetilde{B} - A^{-1} B A^{-1}\widetilde{A} )\\
    $
    \newline
    Below is a translation of the shorthand I use to represent matrices dependent on ``excitation" with respect to the reference Matrix and sums of matrices. Above this line I have represented these excitation matrices with a subscript ``I" but from this point on The subscript will be omitted and it is clear that whenever a matrix depends on $P^{T}_I$ and $Q_{I}$ that this is an excitation matrix. The reference matrix is always $A_{0}$ and is always the Hartree Fock Matrix.
    \newline
    \hfill\break
    $
        Y_{5} = TQ \\
        Y_{6} = (P^{T}TQ)^{-1} = \alpha_{I}^{-1}\\
        Y_{7} = \alpha_{I}^{-1} P^{T} \\
        Y_{11} = \widetilde{M}Q \\
        Y_{23} = P^{T}\widetilde{M}Q \\
        Y_{24} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q \\
        Y_{25} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q\alpha_{I}^{-1} \\
        Y_{26} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q\alpha_{I}^{-1}P^{T}\\
    $
    \newline
    So far you will notice that I have not included up or down arrows to specify what spin the matrices are of. This is because we are calculating the derivative of all up or all down spin orbital rotation parameters at a time. If we are finding the up spin derivatives then any term that is down spin will be constant. The following assumes that we are taking up-spin MO rotation parameter derivatives. Of course the down spin expression can be retrieved by swapping the up and down arrows. I have dubbed any expression with lowercase p prefix as a "precursor" to an expression actually used...
    \newline
    \hfill\break
    $
        \dot{C_{I}} = C_{I}*det(A_{I\downarrow})\\
        \ddot{C_{I}} = C_{I}*\hat{O}det(A_{I\downarrow}) \\
        pK1 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (Q\alpha_{I}^{-1}P^{T}) \\
        pK2 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
        pK3 = \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
        pK4 = \sum_{I=1} \dot{C_{I}} det(A_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
        pK5 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T}) \\
    $
    \newline
    Now these p matrices will be used to make various expressions via BLAS commands.
    \newline
    \hfill\break
    $
        K1T = const0^{-1}*pK1.T =const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (Q\alpha_{I}^{-1}P^{T}T) \\
        TK1T = T.K1T = const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (TQ\alpha_{I}^{-1}P^{T}T)\\ \\
        K2AiB = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}A^{-1}\widetilde{B})\\
        TK2AiB = T.K2AiB = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}A^{-1}\widetilde{B})\\
        K2XA =  const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}X\widetilde{A})\\
        TK2XA = T.K2XA = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}X\widetilde{A})\\ \\
        K2T = \frac{const1}{const0^{2}} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
        TK2T = T.K2T =\frac{const1}{const0^{2}} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}T) \\
        MK2T = \frac{const0}{const1} Y_{4}.K2T= const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (\widetilde{M}Q\alpha_{I}^{-1}P^{T}T)\\ \\
        K3T = const0^{-1}  \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
        TK3T = T.K3T  = const0^{-1}  \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}T)\\ \\
        K4T = \sum_{I=1} \dot{C_{I}} det(A_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
        TK4T = T.K4T = \sum_{I=1} \dot{C_{I}} det(A_{I}) (TQ\alpha_{I}^{-1}P^{T}T) \\ \\
        K5T =  const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T} T)  \\
        TK5T = T.K5T  = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (T Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T} T)  \\
    $
    \newline
    Now with all these matrices and constants the expressions of dhpsioverpsi and dlogpsi can be created.
    In addition I will be using a special generalization of the kinetic operator which I will denote as O. Our Slater matrix with the special O operator applied to each element will be called B_bar

    $
    ``Bbar"_{i,j} =\nabla^2 \phi_{j}(r_i) + \frac{\nabla_{i}J}{J} \cdot \nabla \phi_{j}(r_{i})  + \frac{\nabla^2_i J}{J} \phi_{j}(r_{i})
    $
    */
    ValueMatrix_t Table;
    ValueMatrix_t Bbar;
    ValueMatrix_t Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y11,Y23,Y24,Y25,Y26;
    ValueMatrix_t pK1,K1T,TK1T, pK2,K2AiB,TK2AiB,K2XA,TK2XA,K2T,TK2T,MK2T, pK3,K3T,TK3T, pK4,K4T,TK4T, pK5,K5T,TK5T;

  };
}
#endif
