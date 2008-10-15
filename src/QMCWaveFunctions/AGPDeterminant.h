//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file AGPDeterminant.h
 * @brief Declaration of AGPDeterminant for pairing orbitals.
 */
#ifndef QMCPLUSPLUS_AGP_DIRACDETERMINANT_H
#define QMCPLUSPLUS_AGP_DIRACDETERMINANT_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "QMCWaveFunctions/BasisSetBase.h"

namespace qmcplusplus {

  class AGPDeterminant: public OrbitalBase 
  {

  public:

    ///define BasisSetType with RealType
    typedef BasisSetBase<RealType> BasisSetType;
    typedef BasisSetType::IndexVector_t IndexVector_t;
    typedef BasisSetType::ValueVector_t ValueVector_t;
    typedef BasisSetType::ValueMatrix_t ValueMatrix_t;
    typedef BasisSetType::GradVector_t  GradVector_t;
    typedef BasisSetType::GradMatrix_t  GradMatrix_t;

    BasisSetType* GeminalBasis;

    /** constructor
     *@param spos the single-particle orbital set
     *@param first index of the first particle
     */
    AGPDeterminant(BasisSetType* bs=0);

    ///default destructor
    ~AGPDeterminant();
  
    void checkInVariables(opt_variables_type& active);
    void checkOutVariables(const opt_variables_type& active);
    void resetParameters(const opt_variables_type& active);
    void reportStatus(ostream& os);

    void resetTargetParticleSet(ParticleSet& P);

    ///reset the size: with the number of particles and number of orbtials
    void resize(int nup, int ndown);

    ValueType registerData(ParticleSet& P, PooledData<RealType>& buf);

    ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);

    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

    /** dump the inverse to the buffer
     */
    inline void dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      buf.add(psiM.begin(),psiM.end());
    }

    /** copy the inverse from the buffer
     */
    inline void dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      buf.get(psiM.begin(),psiM.end());
    }

    /** return the ratio only for the  iat-th partcle move
     * @param P current configuration
     * @param iat the particle thas is being moved
     */
    ValueType ratio(ParticleSet& P, int iat);

    /** return the ratio
     * @param P current configuration
     * @param iat particle whose position is moved
     * @param dG differential Gradients
     * @param dL differential Laplacians
     *
     * Data member *_temp contain the data assuming that the move is accepted
     * and are used to evaluate differential Gradients and Laplacians.
     */
    ValueType ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL);

    void ratioUp(ParticleSet& P, int iat);

    void ratioDown(ParticleSet& P, int iat);

    ValueType logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL) {
      ValueType r=ratio(P,iat,dG,dL);
      return evaluateLogAndPhase(r,PhaseValue);
    }


    /** move was accepted, update the real container
     */
    void acceptMove(ParticleSet& P, int iat);

    /** move was rejected. copy the real container to the temporary to move on
     */
    void restore(int iat);

    
    void update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat);


    ValueType evaluateLog(ParticleSet& P, PooledData<RealType>& buf);


    void resizeByWalkers(int nwalkers);

    ///evaluate log of determinant for a particle set: should not be called 
    ValueType
    evaluateLog(ParticleSet& P, 
	        ParticleSet::ParticleGradient_t& G, 
	        ParticleSet::ParticleLaplacian_t& L);
    //{
    //  ValueType psi=evaluate(P,G,L);
    //  return LogValue = evaluateLogAndPhase(psi,PhaseValue);
    //}

    /** Calculate the value of the Dirac determinant for particles
     *@param P input configuration containing N particles
     *@param G a vector containing N gradients
     *@param L a vector containing N laplacians
     *@return the value of the determinant
     *
     *\f$ (first,first+nel). \f$  Add the gradient and laplacian 
     *contribution of the determinant to G(radient) and L(aplacian)
     *for local energy calculations.
     */ 
    ValueType
    evaluate(ParticleSet& P, 
	     ParticleSet::ParticleGradient_t& G, 
	     ParticleSet::ParticleLaplacian_t& L);

    OrbitalBasePtr makeClone(ParticleSet& tqp) const;

    ///Total number of particles
    int NumPtcls;
    ///number of major spins
    int Nup;
    ///number of minor spins
    int Ndown;
    ///size of the basis set
    int BasisSize;

    ///index of the particle (or row) 
    int WorkingIndex;      

    /////Current determinant value
    //ValueType CurrentDet;

    ///coefficient of the up/down block
    ValueMatrix_t Lambda;

    ///coefficient of the major block
    ValueMatrix_t LambdaUP;

    /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
    ValueMatrix_t psiM, psiM_temp;


    /**  Transient data for gradient and laplacian evaluation
     *
     * \f$phiD(j,k) = \sum_{j^{'}} \lambda_{j^{'},j} \phi_k(r_j) \f$
     * j runs over the particle index index
     */
    ValueMatrix_t phiT;

    /// temporary container for testing
    ValueMatrix_t psiMinv;
    /// store gradients
    GradMatrix_t dY;
    /// store laplacians
    ValueMatrix_t d2Y;
    /// temporary determinant-related matrix for gradients
    GradMatrix_t dpsiU, dpsiD;
    /// temporary determinant-related matrix for laplacians
    ValueMatrix_t d2psiU, d2psiD;

    /// value of single-particle orbital for particle-by-particle update
    /** temporary vector for a particle-by-particle move
     *
     * phiTv = Lambda Y(iat)
     */
    ValueVector_t phiTv;
    ValueVector_t psiU, psiD;
    GradVector_t  dpsiUv, dpsiDv;
    ValueVector_t d2psiUv, d2psiDv;
    ValueVector_t workV1, workV2;
    ValueVector_t WorkSpace;
    IndexVector_t Pivot;

    ///current ratio
    RealType curRatio;
    ///cummulate ratio for particle-by-particle update
    RealType cumRatio;
    ///address of  dpsiU[0][0]
    BasisSetType::ValueType *FirstAddressOfdVU;
    ///address of FirstAddressOfdVU+OHMMS_DIM*Nup*Nup
    BasisSetType::ValueType *LastAddressOfdVU;
    ///address of  dpsiD[0][0]
    BasisSetType::ValueType *FirstAddressOfdVD;
    ///address of FirstAddressOfdVD+OHMMS_DIM*Ndown*Nup
    BasisSetType::ValueType *LastAddressOfdVD;
    ///address of myG[0][0]
    ValueType *FirstAddressOfG;
    ///address of FirstAddressOfG+OHMMS_DIM*NumPtcls
    ValueType *LastAddressOfG;
    ///address of dY[0][0]
    BasisSetType::ValueType *FirstAddressOfdY;
    ///address of FirstAddressOfdY+NumPtcls*BasisSize
    BasisSetType::ValueType *LastAddressOfdY;

    ParticleSet::ParticleGradient_t myG, myG_temp;
    ParticleSet::ParticleLaplacian_t myL, myL_temp;
    
    void evaluateLogAndStore(ParticleSet& P);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
