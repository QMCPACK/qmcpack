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
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"


#include "Configuration.h"
#include "QMCHamiltonians/RPAPressure.h"
#include "QMCWaveFunctions/Jastrow/RPAJastrow.h"
#include "Optimize/VarList.h"
#include "ParticleBase/ParticleAttribOps.h"




namespace qmcplusplus {

  typedef QMCHamiltonianBase::Return_t Return_t;
  
  void RPAPressure::resetTargetParticleSet(ParticleSet& P) {
    dPsi[0]->resetTargetParticleSet(P);
    pNorm = 1.0/(P.Lattice.DIM*P.Lattice.Volume);
  };
    
  Return_t RPAPressure::evaluate(ParticleSet& P) {
    vector<OrbitalBase*>::iterator dit(dPsi.begin()), dit_end(dPsi.end());
    tValue=0.0;
    Value=0.0;
    dG = 0.0;
    dL = 0.0;
    while(dit != dit_end) {
      tValue += (*dit)-> evaluateLog(P,dG,dL);
      ++dit;
    }
    tValue *= drsdV;
    ZVCorrection =  -1.0 * (0.5*Sum(dL)+Dot(dG,P.G));
    Value = -ZVCorrection*drsdV;
    
    Press=2.0*P.PropertyList[LOCALENERGY]-P.PropertyList[LOCALPOTENTIAL];
    Press*=pNorm;
    Energy = P.PropertyList[LOCALENERGY];
    
    return 0.0;
  }

  Return_t RPAPressure::evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
    return evaluate(P);
  }

  bool RPAPressure::put(xmlNodePtr cur, ParticleSet& P) {
    MyName = "RPAZVZBP";
    RealType tlen=std::pow(0.75/M_PI*P.Lattice.Volume/static_cast<RealType>(P.getTotalNum()),1.0/3.0);
    drsdV= tlen/(3.0* P.Lattice.Volume);
    RPAJastrow* rpajastrow = new RPAJastrow(P);
    rpajastrow->put(cur);
    dPsi.push_back(rpajastrow);

    return true;
  }


  QMCHamiltonianBase* RPAPressure::makeClone(ParticleSet& P, TrialWaveFunction& psi)
  {
    RPAPressure* myClone = new RPAPressure(P);
    
    vector<OrbitalBase*>::iterator dit(dPsi.begin()), dit_end(dPsi.end());
    while(dit != dit_end) {
      myClone->dPsi.push_back((*dit)->makeClone(P)); 
      ++dit;
    }
    return myClone;
  }

}


/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
* $Id: BareKineticEnergy.h 1581 2007-01-04 16:02:14Z jnkim $ 
***************************************************************************/

