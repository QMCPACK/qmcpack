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
    
    
#ifndef QMCPLUSPLUS_NUMERICALPOTENTIAL_H
#define QMCPLUSPLUS_NUMERICALPOTENTIAL_H
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimLinearSpline.h"

namespace qmcplusplus
{

/** Evaluates the NumericalRadial Potential for a set of source and target particles.
 */
struct NumericalRadialPotential: public QMCHamiltonianBase
{
  ///typedef the grid
  typedef OneDimGridBase<RealType> GridType;
  ///typedef radial potential
  typedef OneDimLinearSpline<RealType> RadialPotentialType;
  ///true if the pair
  bool IsPairPotential;
  ///number of centers
  int Centers;
  ///distance table
  DistanceTableData* d_table;
  ///reference to the center particleset
  ParticleSet& sourcePtcl;
  ///function on a grid
  RadialPotentialType* VofR;

  /** constructor for a two-body potential
   */
  NumericalRadialPotential(ParticleSet& center);

  /** constructor for an one-body potential
   */
  NumericalRadialPotential(ParticleSet& center, ParticleSet& visitor);
  ///destructor
  ~NumericalRadialPotential();

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** Do nothing */
  bool put(xmlNodePtr cur);

  bool get(std::ostream& os) const;

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

};
}
#endif


