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
/** @file DiracDeterminantBase.h
 * @brief declaration of DiracDeterminantBase
 *
 * This is an abstract base class of DiractDeterminant.
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTBASE_H
#define QMCPLUSPLUS_DIRACDETERMINANTBASE_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus {

  class ParticleSet;
  class WalkerSetRef;

  /** Base class to handle a DiracDeterminantBase.
   */
  class DiracDeterminantBase: public OrbitalBase {
  public:
    typedef Matrix<ValueType> Determinant_t;
    typedef Matrix<GradType>  Gradient_t;
    typedef Matrix<ValueType> Laplacian_t;
    typedef std::vector<ValueType> TempValue_t;
    typedef std::vector<GradType>  TempGradient_t;
    typedef std::vector<ValueType> TempLaplacian_t;

    /** constructor
     *@param spos the single-particle orbital set
     *@param first index of the first particle
     */
    DiracDeterminantBase(int first=0);

    ///default destructor
    virtual ~DiracDeterminantBase();
  

    DiracDeterminantBase& operator=(const DiracDeterminantBase& s);

    /** set the index of the first particle in the determinant and reset the size of the determinant
     *@param first index of first particle
     *@param nel number of particles in the determinant
     */
    inline void set(int first, int nel) {
      FirstIndex = first;
      resize(nel,nel);
    }

    ///reset the size: with the number of particles and number of orbtials
    void resize(int nel, int morb);
    ///resize the data with the number of walkers
    void resizeByWalkers(int nwalkers);
    ///return the number of rows (or the number of electrons)
    inline int rows() const { return NumPtcls;}
    ///return the number of coloumns  (or the number of orbitals)
    inline int cols() const { return NumOrbitals;}

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

    ValueType ratio(ParticleSet& P, int iat);
    ValueType ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL);
    ValueType logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL);

    void update(ParticleSet& P, int iat);
    void restore(int iat);
    void update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat);

    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf);

    ValueType
    evaluateLog(ParticleSet& P, 
	        ParticleSet::ParticleGradient_t& G, 
	        ParticleSet::ParticleLaplacian_t& L) {
      std::cerr << "DiracDeterminantBase::evaluateLog should never be called directly" << std::endl;
      return 0.0;
    }


    ValueType
    evaluate(ParticleSet& P, 
	     ParticleSet::ParticleGradient_t& G, 
	     ParticleSet::ParticleLaplacian_t& L);

    virtual void evaluateSingle(ParticleSet& P, int iat)=0;
    virtual void evaluateSingleAll(ParticleSet& P, int iat)=0;
    virtual void evaluateAll(ParticleSet& P)=0;

    ///The number of particles
    int NP;

    ///Number of orbitals
    int NumOrbitals;

    ///Number of particles
    int NumPtcls;

    ///index of the first particle with respect to the particle set
    int FirstIndex;

    ///index of the last particle with respect to the particle set
    int LastIndex;

    ///index of the particle (or row) 
    int WorkingIndex;      

    ///Current determinant value
    ValueType CurrentDet;

    /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
    Determinant_t psiM, psiM_temp;

    /// temporary container for testing
    Determinant_t psiMinv;

    /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$
    Gradient_t    dpsiM, dpsiM_temp;

    /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$
    Laplacian_t   d2psiM, d2psiM_temp;

    /// value of single-particle orbital for particle-by-particle update
    TempValue_t psiV;
    TempGradient_t dpsiV;
    TempLaplacian_t d2psiV;
    TempValue_t workV1, workV2;

    ///storages to process many walkers once
    //vector<Determinant_t> psiM_v; 
    //vector<Gradient_t>    dpsiM_v; 
    //vector<Laplacian_t>   d2psiM_v; 

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
