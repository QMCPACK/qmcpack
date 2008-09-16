//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file DiracDeterminantBaseBase.h
 * @brief Declaration of DiracDeterminantBase with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTWITHBASE_H
#define QMCPLUSPLUS_DIRACDETERMINANTWITHBASE_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"

namespace qmcplusplus {

  class DiracDeterminantBase: public OrbitalBase {
  public:

    typedef SPOSetBase::IndexVector_t IndexVector_t;
    typedef SPOSetBase::ValueVector_t ValueVector_t;
    typedef SPOSetBase::ValueMatrix_t ValueMatrix_t;
    typedef SPOSetBase::GradVector_t  GradVector_t;
    typedef SPOSetBase::GradMatrix_t  GradMatrix_t;

    /** constructor
     *@param spos the single-particle orbital set
     *@param first index of the first particle
     */
    DiracDeterminantBase(SPOSetBasePtr const &spos, int first=0);

    ///default destructor
    ~DiracDeterminantBase();
  
    /**copy constructor
     * @param s existing DiracDeterminantBase
     *
     * This constructor makes a shallow copy of Phi.
     * Other data members are allocated properly.
     */
    DiracDeterminantBase(const DiracDeterminantBase& s);

    DiracDeterminantBase& operator=(const DiracDeterminantBase& s);

    /** return a clone of Phi
     */
    SPOSetBasePtr clonePhi() const;

    inline IndexType rows() const {
      return NumPtcls;
    }

    inline IndexType cols() const {
      return NumOrbitals;
    }

    /** set the index of the first particle in the determinant and reset the size of the determinant
     *@param first index of first particle
     *@param nel number of particles in the determinant
     */
    void set(int first, int nel);

    ///optimizations  are disabled
    inline void checkInVariables(opt_variables_type& active)
    {
      //Phi->checkInVariables(active); 
    }

    inline void checkOutVariables(const opt_variables_type& active)
    {
      //Phi->checkOutVariables(active); 
    }
    
    void resetParameters(const opt_variables_type& active) 
    { 
      //Phi->resetParameters(active); 
    }

    inline void reportStatus(ostream& os)
    {
    }
    void resetTargetParticleSet(ParticleSet& P) { 
      Phi->resetTargetParticleSet(P);
    }

    ///reset the size: with the number of particles and number of orbtials
    void resize(int nel, int morb);

    ValueType registerData(ParticleSet& P, PooledData<RealType>& buf);

    ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);

    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

    /** dump the inverse to the buffer
     */
    void dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf);

    /** copy the inverse from the buffer
     */
    void dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

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

    ValueType logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL);

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


    ///evaluate log of determinant for a particle set: should not be called 
    ValueType
    evaluateLog(ParticleSet& P, 
	        ParticleSet::ParticleGradient_t& G, 
	        ParticleSet::ParticleLaplacian_t& L) ;


    ValueType
    evaluate(ParticleSet& P, 
	     ParticleSet::ParticleGradient_t& G, 
	     ParticleSet::ParticleLaplacian_t& L);

    OrbitalBasePtr makeClone(ParticleSet& tqp) const;

    bool UseRatioOnly;
    ///The number of particles
    int NP;

    ///number of single-particle orbitals which belong to this Dirac determinant
    int NumOrbitals;
    ///number of particles which belong to this Dirac determinant
    int NumPtcls;

    ///index of the first particle with respect to the particle set
    int FirstIndex;

    ///index of the last particle with respect to the particle set
    int LastIndex;

    ///a set of single-particle orbitals used to fill in the  values of the matrix 
    SPOSetBasePtr Phi;

    ///index of the particle (or row) 
    int WorkingIndex;      

    ///Current determinant value
    ValueType CurrentDet;
    ///the sign of determinant
    RealType DetSign;
    /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
    ValueMatrix_t psiM, psiM_temp;

    /// temporary container for testing
    ValueMatrix_t psiMinv;

    /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$
    GradMatrix_t  dpsiM, dpsiM_temp;

    /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$
    ValueMatrix_t d2psiM, d2psiM_temp;

    /// value of single-particle orbital for particle-by-particle update
    ValueVector_t psiV;
    GradVector_t dpsiV;
    ValueVector_t d2psiV;
    ValueVector_t workV1, workV2;

    Vector<ValueType> WorkSpace;
    Vector<IndexType> Pivot;

    ValueType curRatio,cumRatio;
    ValueType *FirstAddressOfG;
    ValueType *LastAddressOfG;
    ValueType *FirstAddressOfdV;
    ValueType *LastAddressOfdV;

    ParticleSet::ParticleGradient_t myG, myG_temp;
    ParticleSet::ParticleLaplacian_t myL, myL_temp;
  };


  
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
