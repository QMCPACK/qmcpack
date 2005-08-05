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
#ifndef OHMMS_QMC_HARMONICPOTENTIAL_H
#define OHMMS_QMC_HARMONICPOTENTIAL_H
#include <algo.h>
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace ohmmsqmc {

  /** Evaluates the Harmonic Potential for a set of source and target particles.
   *
   * \f[ H = \sum_I \frac{1}{2}\omega(I)^2 r^2 \f] 
   * where \f$ \omega(I) \f$ is the frequency of oscillation 
   * around the \f$ Ith \f$ center.
   */
  struct HarmonicPotential: public QMCHamiltonianBase {
    ///number of centers
    int Centers;
    ///container for \f$ 0.5\omega^2 \f$ for each center
    vector<RealType> Omega;
    ///distance table
    DistanceTableData* d_table;

    ///constructor
    HarmonicPotential(RealType omega=1.0): Omega(1,omega){ }

    HarmonicPotential(ParticleSet& center, ParticleSet& visitor) { 
      d_table = DistanceTable::getTable(DistanceTable::add(center,visitor));
      int charge = center.Species.addAttribute("charge");
      Centers = center.getTotalNum();
      Omega.resize(Centers);

      RealType C = 0.5;
      for(int iat=0; iat<Centers;iat++) {
        RealType omega = center.Species(charge,center.GroupID[iat]);
        Omega[iat] = C*omega*omega;
      }
    }

    ///destructor
    ~HarmonicPotential() { }

    inline Return_t 
    evaluate(ParticleSet& P) {
      Value=0.0;
      for(int iat=0; iat<Centers; iat++) {
        RealType e = 0.0;
        for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++) {
          e += d_table->r(0,nn)*d_table->r(0,nn);
        }
        Value += Omega[iat]*e;
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

