#ifndef OHMMS_QMC_POLARISATIONPOTENTIAL_H
#define OHMMS_QMC_POLARISATIONPOTENTIAL_H
#include <algo.h>
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace ohmmsqmc {


  struct PolarizationPotential: public QMCHamiltonianBase {

    RealType Efield;
    PolarizationPotential(double field):Efield(field){ }

    ~PolarizationPotential() { }

    inline ValueType 
    evaluate(ParticleSet& P) {
      RealType sum = 0.0;
      for(int i=0; i < P.getTotalNum(); i++) sum += P.R[i][2];
      return (Efield * sum);
    }

    inline ValueType 
    evaluate(ParticleSet& P, RealType& x){
      RealType sum = 0.0;
      for(int i=0; i < P.getTotalNum(); i++) sum += P.R[i][2];
      x = Efield * sum;
      return x;
    }

#ifdef USE_FASTWALKER
    inline void 
    evaluate(WalkerSetRef& W, ValueVectorType& LE) {
    }
#else
    inline void 
    evaluate(WalkerSetRef& W, ValueVectorType& LE) {
    }
#endif
  };
}
#endif

/************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id: PolarizationPotential.h,v 1.1.1.1 2004/08/24 19:21:11 jnkim
 * Exp $
************************************************************************/

