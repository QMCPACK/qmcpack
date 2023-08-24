//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SHO_BASIS_BUILDERT_H
#define QMCPLUSPLUS_SHO_BASIS_BUILDERT_H

#include "QMCWaveFunctions/HarmonicOscillator/SHOSetT.h"
#include "QMCWaveFunctions/SPOSetBuilderT.h"
#include "QMCWaveFunctions/SPOSetInfo.h"

namespace qmcplusplus
{
template<class T>
class SHOSetBuilderT : public SPOSetBuilderT<T>
{
public:
  using RealType  = typename SPOSetT<T>::RealType;
  using PosType   = typename SPOSetT<T>::PosType;
  using indices_t = typename SPOSetBuilderT<T>::indices_t;

  ParticleSet& Ps;

  RealType length;
  RealType mass;
  RealType energy;
  PosType center;

  int nstates;
  int nmax;
  TinyVector<int, QMCTraits::DIM> ind_dims;

  SPOSetInfoSimple<SHOState> basis_states;

  //construction/destruction
  SHOSetBuilderT(ParticleSet& P, Communicate* comm);

  ~SHOSetBuilderT() override;

  //reset parameters
  void reset();

  //SPOSetBuilder interface
  std::unique_ptr<SPOSetT<T>> createSPOSetFromXML(xmlNodePtr cur) override;

  std::unique_ptr<SPOSetT<T>> createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input) override;

  //local functions
  void update_basis_states(int smax);
  void report(const std::string& pad = "");
};

} // namespace qmcplusplus
#endif
