//////////////////////////////////////////////////////////////////
// (c) Copyright 2004  by Dyutiman Das
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Loomis Laboratory of Physics
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: ddas@uiuc.edu
//   Tel:    217-333-3087 (Office) 
//
// Supported by 
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_QMC_EFFMKINETICENERGY_H
#define OHMMS_QMC_EFFMKINETICENERGY_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "Numerics/MatGrid1D.h"
#include "Numerics/Spline3D/Grid1D.h"

namespace ohmmsqmc {

  struct EffMKineticEnergy: public QMCHamiltonianBase {

    RealType M;
    MatGrid1D* inveffm_grid;

    /// The Constructor
    EffMKineticEnergy(const Grid1D& aGrid1D,
		      const std::vector<int>& intvals,
		      const std::vector<int>& priority,
		      const std::vector<double>& inveffm,
		      RealType m=1.0): M(0.5/m){

      /// assign the grid and initialise it 
      inveffm_grid = new MatGrid1D(intvals.size()-1);
      inveffm_grid->init(aGrid1D,intvals,priority,inveffm);

    }

    /// The destructor
    ~EffMKineticEnergy() { if(inveffm_grid) delete inveffm_grid; }

    ValueType 
    evaluate(ParticleSet& P) {
      RealType ke = 0.0;
      double mterm = 0.7322;   /// gradmeff term // ????? CHANGE later
			    /// !!! ddas

      for(int i=0; i<P.getTotalNum(); i++) {
	double z = P.R[i][2];
	RealType M = -0.5*inveffm_grid->prop(z);
	ke += M*(dot(P.G(i),P.G(i)) + P.L(i));
	if(z >= 1700.748 && z <= 1719.6452) ke += mterm*P.G(i)[2];
      }
      return ke;
    }

    ValueType
    evaluate(ParticleSet& P, RealType& x) {
      return x=evaluate(P);
    }


#ifdef WALKERSTORAGE_COLUMNMAJOR
    void evaluate(WalkerSetRef& P, ValueVectorType& LE) {
      ValueVectorType KE(P.walkers());
      for(int iat=0; iat< P.particles(); iat++) {
	for(int iw=0; iw< P.walkers(); iw++) {
	  KE[iw] -= dot(P.G(iw,iat),P.G(iw,iat)) + P.L(iw,iat);
	}
      }
      LE -= M*KE;
    }
#else
    void evaluate(WalkerSetRef& P, ValueVectorType& LE) {
      for(int iw=0; iw< P.walkers(); iw++) {
	RealType ke = 0.0;
	for(int iat=0; iat< P.particles(); iat++) {
	  ke -= dot(P.G(iw,iat),P.G(iw,iat)) + P.L(iw,iat);
	}
	LE[iw] -= M*ke;
      }
    }
#endif
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

