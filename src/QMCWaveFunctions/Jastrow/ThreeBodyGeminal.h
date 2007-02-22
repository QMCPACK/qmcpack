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
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Optimize/VarList.h"
#include "QMCWaveFunctions/BasisSetBase.h"

namespace qmcplusplus {

  /** @ingroup OrbitalComponent
   * @brief ThreeBodyGeminal functions
   */ 
  class ThreeBodyGeminal: public OrbitalBase {

  public:

    typedef BasisSetBase BasisSetType;

    ///constructor
    ThreeBodyGeminal(ParticleSet& ions, ParticleSet& els);

    ~ThreeBodyGeminal();

    ///reset the value of all the Two-Body Jastrow functions
    void resetParameters(OptimizableSetType& optVariables);
    //evaluate the distance table with els
    void resetTargetParticleSet(ParticleSet& P);

    ValueType evaluateLog(ParticleSet& P,
		          ParticleSet::ParticleGradient_t& G, 
		          ParticleSet::ParticleLaplacian_t& L);

    ValueType evaluate(ParticleSet& P,
		       ParticleSet::ParticleGradient_t& G, 
		       ParticleSet::ParticleLaplacian_t& L) {
      return std::exp(evaluateLog(P,G,L));
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

    void acceptMove(ParticleSet& P, int iat);

    inline void update(ParticleSet& P, 		
		       ParticleSet::ParticleGradient_t& dG, 
		       ParticleSet::ParticleLaplacian_t& dL,
		       int iat);

    ValueType registerData(ParticleSet& P, PooledData<RealType>& buf);

    ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf);
    
    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf);

    void setBasisSet(BasisSetType* abasis) { GeminalBasis=abasis;}

    bool put(xmlNodePtr cur, OptimizableSetType& varlist);

    void addOptimizables(OptimizableSetType& varlist);

  private:

    ///reference to the center
    const ParticleSet& CenterRef;
    ///distance table
    const DistanceTableData* d_table;
    ///size of the localized basis set
    int BasisSize;
    ///number of particles
    int NumPtcls;
    ///offset of the index
    int IndexOffSet;
    /** temporary value for update */
    ValueType diffVal;
    ///root name for Lambda compoenents
    string ID_Lambda;
    /** Y(iat,ibasis) value of the iat-th ortbial, the basis index ibasis
     */
    Matrix<ValueType> Y;
    /** dY(iat,ibasis) value of the iat-th ortbial, the basis index ibasis
     */
    Matrix<GradType>  dY;
    /** d2Y(iat,ibasis) value of the iat-th ortbial, the basis index ibasis
     */
    Matrix<ValueType> d2Y;
    /** V(i,j) = Lambda(k,kk) U(i,kk)
     */
    Matrix<RealType> V;

    /** Symmetric matrix connecting Geminal Basis functions */
    Matrix<RealType> Lambda;

    /** Uk[i] = \sum_j dot(U[i],V[j]) */
    Vector<RealType> Uk;

    /** Gradient for update mode */
    Matrix<GradType> dUk;

    /** Laplacian for update mode */
    Matrix<RealType> d2Uk;

    /** temporary Laplacin for update */
    Vector<ValueType> curLap, tLap;
    /** temporary Gradient for update */
    Vector<GradType> curGrad, tGrad;
    /** tempory Lambda*newY for update */
    Vector<ValueType> curV;
    /** tempory Lambda*(newY-Y(iat)) for update */
    Vector<ValueType> delV;
    /** tempory Lambda*(newY-Y(iat)) for update */
    Vector<ValueType> curVal;

    ValueType *FirstAddressOfdY;
    ValueType *LastAddressOfdY;
    ValueType *FirstAddressOfgU;
    ValueType *LastAddressOfgU;

    /** Geminal basis function */
    BasisSetType *GeminalBasis;

    /** evaluateLog and store data for particle-by-particle update */
    void evaluateLogAndStore(ParticleSet& P);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

