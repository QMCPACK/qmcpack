//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_THREEBODY_GEMINAL_H
#define QMCPLUSPLUS_THREEBODY_GEMINAL_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.h"

namespace qmcplusplus {

  /** @ingroup OrbitalComponent
   * @brief ThreeBodyGeminal functions
   */ 
  class ThreeBodyGeminal: public OrbitalBase {

  public:

    typedef GTOMolecularOrbitals::BasisSetType BasisSetType;

    ///constructor
    ThreeBodyGeminal(ParticleSet& ions, ParticleSet& els);

    ~ThreeBodyGeminal();

    ///reset the value of all the Two-Body Jastrow functions
    void reset();

    //evaluate the distance table with els
    void resetTargetParticleSet(ParticleSet& P) {
      d_table = DistanceTable::getTable(DistanceTable::add(d_table->origin(),P));
    }

    ValueType evaluateLog(ParticleSet& P,
		          ParticleSet::ParticleGradient_t& G, 
		          ParticleSet::ParticleLaplacian_t& L);

    ValueType evaluate(ParticleSet& P,
		       ParticleSet::ParticleGradient_t& G, 
		       ParticleSet::ParticleLaplacian_t& L) {
      return exp(evaluateLog(P,G,L));
    }

    ValueType ratio(ParticleSet& P, int iat);

    /** later merge the loop */
    ValueType ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL);

    /** later merge the loop */
    ValueType logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL);

    void restore(int iat);

    void update(ParticleSet& P, int iat);

    inline void update(ParticleSet& P, 		
		       ParticleSet::ParticleGradient_t& dG, 
		       ParticleSet::ParticleLaplacian_t& dL,
		       int iat);

    ValueType registerData(ParticleSet& P, PooledData<RealType>& buf);

    ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf);
    
    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf);

    void setBasisSet(BasisSetType* abasis) { GeminalBasis=abasis;}

    bool put(xmlNodePtr cur, VarRegistry<RealType>& varlist);

  private:

    const DistanceTableData* d_table;
    int BasisSize;
    int NumPtcls;
    ParticleSet& CenterRef;
    /** U(i,k)  k-th element of the i-th particle 
     */
    Matrix<RealType> U;
    /** V(i,j) = Lambda(k,kk) U(i,kk)
     */
    Matrix<RealType> V;

    /** Symmetric matrix connecting Geminal Basis functions */
    Matrix<RealType> Lambda;

    /** Uk[i] = \sum_j dot(U[i],V[j]) */
    vector<RealType> Uk;

    /** Gradient for update mode */
    vector<GradType> dUk;

    /** Laplacian for update mode */
    vector<RealType> d2Uk;

    /** temporary value for update */
    ValueType curVal;
    /** temporary Laplacin for update */
    ValueType curLap;
    /** temporary Gradient for update */
    GradType curGrad;

    ValueType *FirstAddressOfgU;
    ValueType *LastAddressOfgU;

    /** Geminal basis function */
    BasisSetType *GeminalBasis;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

