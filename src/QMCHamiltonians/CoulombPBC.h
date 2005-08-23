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
#include "LongRange/LPQHIBasis.h"

namespace ohmmsqmc {

  /** @ingroup hamiltonian
   *\brief Calculates the AA Coulomb potential using PBCs
   */

  struct CoulombPBC: public QMCHamiltonianBase {

    bool FirstTime;
    ParticleSet& PtclRef;
    LRCoulombAA<LPQHIBasis>* AA;
    
    CoulombPBC(ParticleSet& ref): FirstTime(true), PtclRef(ref), 
				  AA(0) {}
    
    ~CoulombPBC() { 
      if(AA) delete AA;
    }

    inline Return_t evaluate(ParticleSet& P) {  

      if(FirstTime || PtclRef.tag() == P.tag()) {
        //d_ii->evaluate(PtclRef);
        Value = 0.0;

	if(FirstTime)
	  AA = new LRCoulombAA<LPQHIBasis>(PtclRef);

        Value = AA->evalTotal();

	cout << "Finished EvalTotal " << endl;

	cout << Value << endl;
	exit(0);

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

