//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_FREE_ORBITAL_BUILDERT_H
#define QMCPLUSPLUS_FREE_ORBITAL_BUILDERT_H

#include "QMCWaveFunctions/SPOSetBuilderT.h"

namespace qmcplusplus
{
template<typename T>
class FreeOrbitalBuilderT : public SPOSetBuilderT<T>
{
public:
  using RealType = typename SPOSetBuilderT<T>::RealType;
  using PosType  = typename SPOSetBuilderT<T>::PosType;

  FreeOrbitalBuilderT(ParticleSetT<T>& els, Communicate* comm, xmlNodePtr cur);
  ~FreeOrbitalBuilderT() {}

  std::unique_ptr<SPOSetT<T>> createSPOSetFromXML(xmlNodePtr cur) override;

private:
  ParticleSetT<T>& targetPtcl;
  bool in_list(const int j, const std::vector<int> l);
};
} // namespace qmcplusplus
#endif
