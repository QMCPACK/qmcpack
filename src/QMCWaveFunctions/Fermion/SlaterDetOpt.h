///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file src/QMCWaveFunctions/Fermion/SlaterDetOpt.h
///
/// \brief   A class for a Slater determinant with optimizable orbitals.
///
/// \author  Eric Neuscamman
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SLATERDETOPT_H
#define QMCPLUSPLUS_SLATERDETOPT_H

#include <QMCWaveFunctions/OrbitalBase.h>
#include <QMCWaveFunctions/FermionBase.h>
#include <QMCWaveFunctions/SPOSetBase.h>

namespace qmcplusplus {

class TrialWaveFunction;

class LCOrbitalSetOptTrialFunc;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  A class for a Slater determinant with optimizable orbitals.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
class SlaterDetOpt: public OrbitalBase, public FermionBase {

  // private data members
  private:

    /// \brief   pointer to the set of optimizable single particle orbitals
    SPOSetBase * m_spo;

    /// \brief   whether this is for up or down spins (0 for up, 1 for down)
    int m_up_or_down;

    /// \brief   index of the determinant's first electron
    int m_first;

    /// \brief   index one past the determinant's last electron
    int m_last;

    /// \brief   number of electrons
    int m_nel;

    /// \brief   total number of molecular orbitals (i.e. linear combinations) in the optimizable set, including those not occupied in this determinant
    int m_nmo;

    // Ratio of new to old values of the wave function, after a particle move.
    ValueType curRatio;

    /// \brief   matrix of orbital values for orbitals in this determinant
    SPOSetBase::ValueMatrix_t m_orb_val_mat;

    /// \brief   inverse of the orbital value matrix
    SPOSetBase::ValueMatrix_t m_orb_inv_mat;

    /// \brief   matrix of orbital gradients for orbitals in this determinant (each element is a length 3 tiny vector)
    SPOSetBase::GradMatrix_t  m_orb_der_mat;

    /// \brief   matrix of x,y,z summed orbital laplacians for orbitals in this determinant (each element is the sum of the d2dx2, d2dy2, and d2dz2 derivatives for the orbital)
    SPOSetBase::ValueMatrix_t m_orb_lap_mat;

    /// \brief   matrix of orbital values for all molecular orbitals
    SPOSetBase::ValueMatrix_t m_orb_val_mat_all;

    /// \brief   matrix of orbital gradients for all molecular orbitals (each element is a length 3 tiny vector)
    SPOSetBase::GradMatrix_t  m_orb_der_mat_all;

    /// \brief   matrix of x,y,z summed orbital laplacians for all molecular orbitals (each element is the sum of the d2dx2, d2dy2, and d2dz2 derivatives for the orbital)
    SPOSetBase::ValueMatrix_t m_orb_lap_mat_all;

    /// \brief   vector of orbital values for orbitals in this determinant for a particular particle
    SPOSetBase::ValueVector_t m_orb_val_vec;

    /// \brief   vector of orbital gradients for orbitals in this determinant for a particular particle
    SPOSetBase::GradVector_t  m_orb_der_vec;

    /// \brief   vector of x,y,z summed orbital laplacians for orbitals in this determinant for a particular particle
    SPOSetBase::ValueVector_t m_orb_lap_vec;

    /// \brief   matrix to hold partial derivatives of the log of determinant with respect to molecular orbital values
    SPOSetBase::ValueMatrix_t m_dp0;

    /// \brief   matrix to hold partial derivatives of the local energy with respect to molecular orbital values
    SPOSetBase::ValueMatrix_t m_dh0;

    /// \brief   matrix to hold partial derivatives of the local energy with respect to molecular orbital gradients
    SPOSetBase::ValueMatrix_t m_dh1;

    /// \brief   matrix to hold partial derivatives of the local energy with respect to molecular orbital laplacians
    SPOSetBase::ValueMatrix_t m_dh2;

    /// \brief   workspace
    std::vector<RealType> m_work;

    /// \brief   pivot workspace
    std::vector<int> m_pivot;

  // private member functions
  private:

    OrbitalBase::RealType evaluate_matrices_from_scratch(ParticleSet& P, const bool all);

  // public type definitions
  public:

    typedef OrbitalSetTraits<ValueType>::IndexVector_t IndexVector_t;
    typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;
    typedef OrbitalSetTraits<ValueType>::GradVector_t  GradVector_t;
    typedef OrbitalSetTraits<ValueType>::HessMatrix_t  HessMatrix_t;
    typedef OrbitalSetTraits<ValueType>::HessType      HessType;
    typedef Array<HessType,3>                          HessArray_t;
    typedef TinyVector<HessType, OHMMS_DIM>            GGGType;
    typedef Vector<GGGType>                            GGGVector_t;
    typedef Matrix<GGGType>                            GGGMatrix_t;
    typedef ParticleSet::Walker_t                      Walker_t;

  // public member functions
  public:

    SlaterDetOpt(ParticleSet & targetPtcl, SPOSetBase * spo_ptr, const int up_or_down);

    ~SlaterDetOpt();

    void add_orbs_to_tf(TrialWaveFunction & twf, const std::string & name);

    void set_spo_optimizable_rotations();

    LCOrbitalSetOptTrialFunc * tfc_ptr(const std::string & calling_func);

    void checkInVariables(opt_variables_type& active);

    void checkOutVariables(const opt_variables_type& active);

    void resetParameters(const opt_variables_type& active);

    void reportStatus(std::ostream& os);

    void resetTargetParticleSet(ParticleSet& P);

    ValueType evaluate(ParticleSet& P ,ParticleSet::ParticleGradient_t& G ,ParticleSet::ParticleLaplacian_t& L);

    RealType evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

    RealType evaluateLog(ParticleSet& P,
                         ParticleSet::ParticleGradient_t& G,
                         ParticleSet::ParticleLaplacian_t& L,
                         PooledData<RealType>& buf,
                         bool fillBuffer);

    GradType evalGrad(ParticleSet& P, int iat);

    ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);

    ValueType ratio(ParticleSet& P, int iat, ParticleSet::ParticleGradient_t& dG,ParticleSet::ParticleLaplacian_t& dL);

    ValueType ratio(ParticleSet& P, int iat);

    void acceptMove(ParticleSet& P, int iat);

    void restore(int iat);

    void update(ParticleSet& P, ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL, int iat);

    RealType evaluateLog(ParticleSet& P,BufferType& buf);

    RealType registerData(ParticleSet& P, BufferType& buf);

    void registerDataForDerivatives(ParticleSet& P, BufferType& buf, int storageType=0);

    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false);

    void copyFromBuffer(ParticleSet& P, BufferType& buf);

    OrbitalBasePtr makeClone(ParticleSet& tqp) const;

    void evaluateDerivatives(ParticleSet& P,
                             const opt_variables_type& optvars,
                             std::vector<RealType>& dlogpsi,
                             std::vector<RealType>& dhpsioverpsi);

    virtual void memoryUsage_DataForDerivatives(ParticleSet& P,long& orbs_only,long& orbs, long& invs, long& dets)
    {
      //Dets[0]->memoryUsage_DataForDerivatives(P,orbs_only,orbs,invs,dets);
      //Dets[1]->memoryUsage_DataForDerivatives(P,orbs_only,orbs,invs,dets);
    }

};

}

#endif
