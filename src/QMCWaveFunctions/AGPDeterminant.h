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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_AGP_DIRACDETERMINANT_H
#define QMCPLUSPLUS_AGP_DIRACDETERMINANT_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
//#include "QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.h"
//#include "QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.h"

namespace qmcplusplus {

  class AGPDeterminant: public OrbitalBase {

  public:

    typedef GridMolecularOrbitals::BasisSetType BasisSetType;
    //typedef STOMolecularOrbitals::BasisSetType BasisSetType;
    //typedef GTOMolecularOrbitals::BasisSetType BasisSetType;
    typedef Matrix<ValueType> Determinant_t;
    typedef Matrix<GradType>  Gradient_t;
    typedef Matrix<ValueType> Laplacian_t;

    BasisSetType* GeminalBasis;

    /** constructor
     *@param spos the single-particle orbital set
     *@param first index of the first particle
     */
    AGPDeterminant(BasisSetType* bs=0);

    ///default destructor
    ~AGPDeterminant();
  
    ///reset the single-particle orbital set
    void reset() { GeminalBasis->reset(); }
   
    void resetTargetParticleSet(ParticleSet& P);

    ///reset the size: with the number of particles and number of orbtials
    void resize(int nup, int ndown);

    ValueType registerData(ParticleSet& P, PooledData<RealType>& buf);

    ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf);

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


    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf);


    void resizeByWalkers(int nwalkers);

    ///evaluate log of determinant for a particle set: should not be called 
    ValueType
    evaluateLog(ParticleSet& P, 
	        ParticleSet::ParticleGradient_t& G, 
	        ParticleSet::ParticleLaplacian_t& L) {
      ValueType psi=evaluate(P,G,L);
      return LogValue = evaluateLogAndPhase(psi,PhaseValue);
    }

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

    bool UseRatioOnly;

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

    ///Current determinant value
    ValueType CurrentDet;

    ///coefficient of the up/down block
    Determinant_t Lambda;

    ///coefficient of the major block
    Determinant_t LambdaUP;

    /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
    Determinant_t psiM, psiM_temp;


    /**  Transient data for gradient and laplacian evaluation
     *
     * \f$phiD(j,k) = \sum_{j^{'}} \lambda_{j^{'},j} \phi_k(r_j) \f$
     * j runs over the particle index index
     */
    Determinant_t phiT;

    /// temporary container for testing
    Determinant_t psiMinv;
    /// store gradients
    Gradient_t dY;
    /// store laplacians
    Laplacian_t d2Y;
    /// temporary determinant-related matrix for gradients
    Gradient_t dpsiU, dpsiD;
    /// temporary determinant-related matrix for laplacians
    Laplacian_t d2psiU, d2psiD;

    /// value of single-particle orbital for particle-by-particle update
    /** temporary vector for a particle-by-particle move
     *
     * phiTv = Lambda Y(iat)
     */
    Vector<ValueType> phiTv;
    Vector<ValueType> psiU, psiD;
    Vector<GradType> dpsiUv, dpsiDv;
    Vector<ValueType> d2psiUv, d2psiDv;
    Vector<ValueType> workV1, workV2;
    Vector<ValueType> WorkSpace;
    Vector<IndexType> Pivot;

    ValueType curRatio,cumRatio;
    ///address of  dpsiU[0][0]
    ValueType *FirstAddressOfdVU;
    ///address of FirstAddressOfdVU+OHMMS_DIM*Nup*Nup
    ValueType *LastAddressOfdVU;
    ///address of  dpsiD[0][0]
    ValueType *FirstAddressOfdVD;
    ///address of FirstAddressOfdVD+OHMMS_DIM*Ndown*Nup
    ValueType *LastAddressOfdVD;
    ///address of myG[0][0]
    ValueType *FirstAddressOfG;
    ///address of FirstAddressOfG+OHMMS_DIM*NumPtcls
    ValueType *LastAddressOfG;
    ///address of dY[0][0]
    ValueType *FirstAddressOfdY;
    ///address of FirstAddressOfdY+NumPtcls*BasisSize
    ValueType *LastAddressOfdY;

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
