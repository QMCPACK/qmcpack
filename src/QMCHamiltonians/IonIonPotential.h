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

  /**
   *\brief Calculates the Ion-Ion potential.
   *
   * \f[ H = \sum_{I<J} \frac{Z(I)Z(J)}{R_{IJ}}, \f]
   * where \f$ Z(I) \f$ is the effective charge of 
   * ion I.
  */

  struct IonIonPotential: public QMCHamiltonianBase {

    ValueType d_sum;
    IonIonPotential(RealType z=1.0): d_sum(0.0) { }

    IonIonPotential(ParticleSet& ref): d_sum(0.0) { 

      IndexType iii = DistanceTable::add(ref);
      DistanceTableData* d_ii = DistanceTable::getTable(iii);
      d_ii->create(1);

      SpeciesSet& tspecies(ref.getSpeciesSet());
      int charge = tspecies.addAttribute("charge");
      int nat = ref.getTotalNum();
      vector<RealType> Z(nat);
      for(int iat=0; iat<nat;iat++) {
        Z[iat] = tspecies(charge,ref.GroupID[iat]);
      }
      d_ii->evaluate(ref);
      d_sum = 0.0;
      for(int iat=0; iat< nat; iat++) {
        RealType esum = 0.0;
	for(int nn=d_ii->M[iat], jat=iat+1; nn<d_ii->M[iat+1]; nn++,jat++) {
          esum += Z[jat]*d_ii->rinv(nn);
	}
        d_sum += esum*Z[iat];
      }
      LOGMSG("Energy of ion-ion interaction " << d_sum)
    }
    
    ~IonIonPotential() { }

    inline ValueType evaluate(ParticleSet& P) {return d_sum;}

    inline ValueType evaluate(ParticleSet& P, RealType& x) {
      //only for consistency
      return x=d_sum;
    }
    
    //inline void 
    //evaluate(WalkerSetRef& P, ValueVectorType& LE) {
    //  LE += d_sum;
    //}
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

