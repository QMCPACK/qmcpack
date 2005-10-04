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
#ifndef OHMMS_QMC_IONIONPOTENTIAL_H
#define OHMMS_QMC_IONIONPOTENTIAL_H
#include <algo.h>
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace ohmmsqmc {

  /** @ingroup hamiltonian
   *\brief Calculates the Ion-Ion potential.
   *
   * \f[ H = \sum_{I<J} \frac{Z(I)Z(J)}{R_{IJ}}, \f]
   * where \f$ Z(I) \f$ is the effective charge of 
   * ion I.
   * @todo IonIonPotential and CoulombPotentialAA should be 
   * merged to one.
   */

  struct IonIonPotential: public QMCHamiltonianBase {

    bool FirstTime;
    DistanceTableData* d_ii;
    ParticleSet& PtclRef;
    vector<RealType> Z;
    
    IonIonPotential(ParticleSet& ref): FirstTime(true), d_ii(0), PtclRef(ref){ 
      
      d_ii = DistanceTable::getTable(DistanceTable::add(ref));

      SpeciesSet& tspecies(ref.getSpeciesSet());
      int charge = tspecies.addAttribute("charge");
      int nat = ref.getTotalNum();
      Z.resize(nat);
      for(int iat=0; iat<nat;iat++) {
        Z[iat] = tspecies(charge,ref.GroupID[iat]);
      }
     
    }
    
    ~IonIonPotential() { }

    inline Return_t evaluate(ParticleSet& P) {  
      //Later it should check if the d_ii is already updated
      if(FirstTime) {
        //d_ii->evaluate(PtclRef);
        Value = 0.0;
        for(int iat=0; iat< Z.size(); iat++) {
          RealType esum = 0.0;
          for(int nn=d_ii->M[iat], jat=iat+1; nn<d_ii->M[iat+1]; nn++,jat++) {
            esum += Z[jat]*d_ii->rinv(nn);
          }
          Value += esum*Z[iat];
        }
        LOGMSG("Energy of ion-ion interaction " << Value)
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

