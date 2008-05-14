//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_RPAPRESSURECORRECTION_H
#define QMCPLUSPLUS_RPAPRESSURECORRECTION_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "LongRange/LRHandlerBase.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"


#include "Configuration.h"

namespace qmcplusplus {

  /** @ingroup hamiltonian
   @brief Evaluate the Bare Pressure
  **/

  struct RPAPressureCorrection: public QMCHamiltonianBase {
    
    typedef CubicBspline<RealType,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> SplineEngineType;
    typedef CubicSplineSingle<RealType,SplineEngineType> FuncType;
    typedef LinearGrid<RealType> GridType;
    typedef LRHandlerBase HandlerType;
    typedef RealType Return_t;
    
    ///Laplacian and Gradient for derivative of wave function with respect to r_s
    PtclOnLatticeTraits::ParticleGradient_t dG;
    PtclOnLatticeTraits::ParticleLaplacian_t dL;
    vector<OrbitalBase*> dPsi;
    Return_t Rs;
    Return_t tValue;
    Return_t drsdV;
    Return_t pNorm;
    Return_t ZVCorrection;
    Return_t drsStar;
    bool ZV;
    bool ZB;

    /** constructor
     *
     * Pressure operators need to be re-evaluated during optimization.
     */
    RPAPressureCorrection(ParticleSet& P): dG(P.G),dL(P.L),Rs(0.0),ZV(false),ZB(false),tValue(0.0),realHandler(0),drsdV(0.0),ZVCorrection(0.0),drsStar(1.0), OwnHandler(true) { 
      UpdateMode.set(OPTIMIZABLE,1);
/*      pNorm = 1.0/(P.Lattice.DIM*P.Lattice.Volume);*/
    };
    
    ///destructor
    ~RPAPressureCorrection() {};

    void resetTargetParticleSet(ParticleSet& P);

    inline Return_t 
    evaluate(ParticleSet& P);

    inline Return_t 
    evaluate(ParticleSet& P, vector<NonLocalData>& Txy) ;
    
    /** implements the virtual function.
     * 
     * Nothing is done but should check the mass
     */
    bool put(xmlNodePtr cur ) {return true;}
    
    bool put(xmlNodePtr cur, ParticleSet& P);

    bool get(std::ostream& os) const {
      os << "ZVZBPressure";
      return true;
    };

    QMCHamiltonianBase* clone(ParticleSet& P, TrialWaveFunction& psi);
    
    private:
//       ShortRangePartAdapter<RealType>* sra;
      HandlerType* realHandler;
      string MyName;
      bool OwnHandler;
      
      /** container to keep track unique instantiations of HandlerType 
       *
       * The key uses jastrow/@name. This object is used to perform the breakup only once
       * for the handler with the same parameters (name). 
      */
      static map<string,HandlerType*> handlerSet;
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: BareKineticEnergy.h 1581 2007-01-04 16:02:14Z jnkim $ 
 ***************************************************************************/

