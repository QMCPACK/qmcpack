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
#ifndef OHMMS_QMC_COULOMBPBC_H
#define OHMMS_QMC_COULOMBPBC_H
#include <algo.h>
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

//Long-range includes.
#include "LongRange/LRCoulombAA.h"
#include "LongRange/LRCoulombAB.h"
#include "LongRange/LPQHIBasis.h"

namespace ohmmsqmc {

  /** @ingroup hamiltonian
   *\brief Calculates the AA Coulomb potential using PBCs
   */

  struct CoulombPBCAA: public QMCHamiltonianBase {

    bool FirstTime;
    //TODO : MAKE THIS A POINTER SO IT CAN BE UPDATED WHEN PARTICLESET CHANGES
    ParticleSet& PtclRef;
    LRCoulombAA<LPQHIBasis>* AA;
    
    CoulombPBCAA(ParticleSet& ref): FirstTime(true), PtclRef(ref), 
				  AA(0) {
      //Create the distancetable here so that it's added to the auto-update
      //list. Pointer is retrieved again in LRCoulombAA.
      /*
      cout << endl << endl;
      LOGMSG("REMOVE THESE AFTER CVS UPDATE");
      cout << endl << endl;

      DistanceTable::add(PtclRef);
      */
    }
    
    ~CoulombPBCAA() { 
      if(AA) delete AA;
    }

    void resetTargetParticleSet(ParticleSet& P) {
      ERRORMSG("CoulombPBCAA::resetTargetParticleSet IS NOT VALID");
      //AA->resetTargetParticleSet(P);
    }

    inline Return_t evaluate(ParticleSet& P) {  

      if(FirstTime || PtclRef.tag() == P.tag()) {
        Value = 0.0;

	if(FirstTime) { //Init the breakup on first call.
	  LOGMSG("Performing long-range breakup for CoulombAA potential");
	  AA = new LRCoulombAA<LPQHIBasis>(PtclRef);
	}

        Value = AA->evalTotal();

	FirstTime = false;
      }
      return Value;
    }

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      return true;
    }
  };


  struct CoulombPBCAB: public QMCHamiltonianBase {

    bool FirstTime;
    ParticleSet &PtclIons, &PtclElns;
    LRCoulombAB<LPQHIBasis>* AB;
    
    CoulombPBCAB(ParticleSet& ions,ParticleSet& elns): FirstTime(true), PtclIons(ions), PtclElns(elns), AB(0) {
      //Create the distancetable here so that it's added to the auto-update
      //list. Pointer is retrieved again in LRCoulombAB.
      //DistanceTable::add(ions,elns);
    }
    
    ~CoulombPBCAB() { 
      if(AB) delete AB;
    }

    inline Return_t evaluate(ParticleSet& P) {  

      if(FirstTime || PtclIons.tag() == P.tag() || PtclElns.tag() == P.tag()) {
        Value = 0.0;

	if(FirstTime) {//Init the breakup on first call.
	  LOGMSG("Performing long-range breakup for CoulombAB potential");
	  AB = new LRCoulombAB<LPQHIBasis>(PtclIons,PtclElns);
	}

        Value = AB->evalTotal();

	FirstTime = false;
      }
      return Value;
    }

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      return true;
    }
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

