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
    
    
#ifndef QMCPLUSPLUS_ONEBODYDENSITYMATRIX_HAMILTONIAN_H
#define QMCPLUSPLUS_ONEBODYDENSITYMATRIX_HAMILTONIAN_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{

class OneBodyDensityMatrix: public QMCHamiltonianBase
{
public:
  OneBodyDensityMatrix(ParticleSet& ions, ParticleSet& elns);

  OneBodyDensityMatrix(ParticleSet& elns);

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist);
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);
  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

private:
  /// maximum distance
  RealType Dmax;
  /// bin size
  RealType Delta;
  /// one of bin size
  RealType DeltaInv;
  ///trial wavefunction
  TrialWaveFunction& Psi;
  ///data
  Vector<RealType> gofr;
  /** resize the internal data
   *
   * The argument list is not completed
   */
  void resize();
};

}
#endif

