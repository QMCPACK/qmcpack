//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
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

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: NumericalRadialPotential.h 1581 2007-01-04 16:02:14Z jnkim $
 ***************************************************************************/

