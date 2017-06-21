//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_IONIONPOTENTIAL_H
#define QMCPLUSPLUS_IONIONPOTENTIAL_H
#if (__GNUC__ == 2)
#include <algo.h>
#else
#include <algorithm>
#endif
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

#ifdef QMC_CUDA
class MCWalkerConfiguration;
#endif

namespace qmcplusplus
{

/** @ingroup hamiltonian
 *\brief Calculates the Ion-Ion potential.
 *
 * \f[ H = \sum_{I<J} \frac{Z(I)Z(J)}{R_{IJ}}, \f]
 * where \f$ Z(I) \f$ is the effective charge of
 * ion I.
 * @todo IonIonPotential and CoulombPotentialAA should be
 * merged to one.
 */

struct IonIonPotential: public QMCHamiltonianBase
{

  bool FirstTime;
  DistanceTableData* d_ii;
  ParticleSet& PtclRef;
  std::vector<RealType> Z;

  IonIonPotential(ParticleSet& ref): FirstTime(true), d_ii(0), PtclRef(ref)
  {
    d_ii = DistanceTable::add(ref);
    SpeciesSet& tspecies(ref.getSpeciesSet());
    int charge = tspecies.addAttribute("charge");
    int nat = ref.getTotalNum();
    Z.resize(nat);
    for(int iat=0; iat<nat; iat++)
    {
      Z[iat] = tspecies(charge,ref.GroupID[iat]);
    }
  }

  ~IonIonPotential() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    //do nothing
  }

  inline Return_t evaluate(ParticleSet& P)
  {
    //Later it should check if the d_ii is already updated
    if(FirstTime)
    {
      Value = 0.0;
      for(int iat=0; iat< Z.size(); iat++)
      {
        RealType esum = 0.0;
        for(int nn=d_ii->M[iat], jat=iat+1; nn<d_ii->M[iat+1]; nn++,jat++)
        {
          esum += Z[jat]*d_ii->rinv(nn);
        }
        Value += esum*Z[iat];
      }
      FirstTime = false;
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
    os << "IonIonPotential: " << PtclRef.getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new IonIonPotential(PtclRef);
  }

  /*@{
   * @brief functions to handle particle-by-particle move
   */
  Return_t registerData(ParticleSet& P, BufferType& buffer)
  {
    return evaluate(P);
  }
  void copyFromBuffer(ParticleSet& P, BufferType& buf) { }
  void copyToBuffer(ParticleSet& P, BufferType& buf) { }
  Return_t evaluatePbyP(ParticleSet& P, int active)
  {
    return Value;
  }
  void acceptMove(int active) { }
  void rejectMove(int active) { }

#ifdef QMC_CUDA
  void addEnergy(MCWalkerConfiguration &W, std::vector<RealType> &LocalEnergy)
  {
    if (FirstTime)
      evaluate(PtclRef);
    std::vector<Walker_t*> &walkers = W.WalkerList;
    for (int iw=0; iw<walkers.size(); iw++)
    {
      walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] = Value;
      LocalEnergy[iw] += Value;
    }
  }
#endif
  /*@}*/
};
}
#endif


