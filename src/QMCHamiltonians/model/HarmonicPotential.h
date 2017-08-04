//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_HARMONICPOTENTIAL_H
#define QMCPLUSPLUS_HARMONICPOTENTIAL_H
#include <algo.h>
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{

/** Evaluates the Harmonic Potential for a set of source and target particles.
 *
 * \f[ H = \sum_I \frac{1}{2}\omega(I)^2 r^2 \f]
 * where \f$ \omega(I) \f$ is the frequency of oscillation
 * around the \f$ Ith \f$ center.
 */
struct HarmonicPotential: public QMCHamiltonianBase
{
  ///number of centers
  int Centers;
  ///container for \f$ 0.5\omega^2 \f$ for each center
  std::vector<RealType> Omega;
  ///distance table
  DistanceTableData* d_table;
  ///reference to the center particleset
  ParticleSet& sourcePtcl;

  HarmonicPotential(ParticleSet& center, ParticleSet& visitor): sourcePtcl(center)
  {
    d_table = DistanceTable::add(center,visitor);
    //int charge = center.Species.addAttribute("charge");
    Centers = center.getTotalNum();
    Omega.resize(Centers,0.5);
    RealType C = 0.5;
    for(int iat=0; iat<Centers; iat++)
    {
      //RealType omega = center.getSpeciesSet(charge,center.GroupID[iat]);
      RealType omega=1.0;
      Omega[iat] = C*omega*omega;
    }
  }

  ///destructor
  ~HarmonicPotential() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    d_table = DistanceTable::add(sourcePtcl,P);
  }

  inline Return_t
  evaluate(ParticleSet& P)
  {
    Value=0.0;
    for(int iat=0; iat<Centers; iat++)
    {
      RealType e = 0.0;
      for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++)
        e += d_table->r(nn)*d_table->r(nn);
      Value += Omega[iat]*e;
    }
    return Value;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }
  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "Harmonic potential: ";
    for(int i=0; i<Centers; i++)
      os << Omega[i] << " ";
    os << std::endl;
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new HarmonicPotential(sourcePtcl,qp);
  }

};
}
#endif


