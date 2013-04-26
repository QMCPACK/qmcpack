#ifndef QMCPLUSPLUS_POLARISATIONPOTENTIAL_H
#define QMCPLUSPLUS_POLARISATIONPOTENTIAL_H
#include <algo.h>
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{


struct PolarizationPotential: public QMCHamiltonianBase
{

  RealType Efield;
  PolarizationPotential(double field):Efield(field) { }

  ~PolarizationPotential() { }

  inline Return_t
  evaluate(ParticleSet& P)
  {
    RealType sum = 0.0;
    for(int i=0; i < P.getTotalNum(); i++)
      sum += P.R[i][2];
    return Value=Efield * sum;
  }

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

};
}
#endif

/************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id: PolarizationPotential.h,v 1.1.1.1 2004/08/24 19:21:11 jnkim
 * Exp $
************************************************************************/

