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
    
    
#ifndef QMCPLUSPLUS_EFFMKINETICENERGY_H
#define QMCPLUSPLUS_EFFMKINETICENERGY_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "Numerics/MatGrid1D.h"
#include "Numerics/Spline3D/Grid1D.h"

namespace qmcplusplus
{

/** Class to evaluate kinectic energy of a heterostructure
 *
 * Spatially varying effective-mass is used to evaluate
 * the kinetic energy term.
 */
struct EffMKineticEnergy: public QMCHamiltonianBase
{

  RealType M;
  MatGrid1D* inveffm_grid;

  /// The Constructor
  EffMKineticEnergy(const Grid1D& aGrid1D,
                    const std::vector<int>& intvals,
                    const std::vector<int>& priority,
                    const std::vector<double>& inveffm,
                    RealType m=1.0): M(0.5/m)
  {
    // assign the grid and initialise it
    inveffm_grid = new MatGrid1D(intvals.size()-1);
    inveffm_grid->init(aGrid1D,intvals,priority,inveffm);
  }

  /// The destructor
  ~EffMKineticEnergy()
  {
    if(inveffm_grid)
      delete inveffm_grid;
  }

  Return_t
  evaluate(ParticleSet& P)
  {
    RealType ke = 0.0;
    double mterm = 0.7322;   // gradmeff term // ????? CHANGE later
    for(int i=0; i<P.getTotalNum(); i++)
    {
      RealType z = P.R[i][2];
      RealType M = -0.5*inveffm_grid->prop(z);
      ke += M*(dot(P.G(i),P.G(i)) + P.L(i));
      if(z >= 1700.748 && z <= 1719.6452)
        ke += mterm*P.G(i)[2];
    }
    return Value=ke;
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
};
}
#endif


