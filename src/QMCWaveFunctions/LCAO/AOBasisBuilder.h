//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_ATOMICORBITALBUILDER_H
#define QMCPLUSPLUS_ATOMICORBITALBUILDER_H


#include "Message/MPIObjectBase.h"
#include "hdf/hdf_archive.h"
#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{
/** atomic basisset builder
   * @tparam COT, CenteredOrbitalType = SoaAtomicBasisSet<RF,SH>
   *
   * Reimplement AtomiSPOSetBuilder.h
   */
template<typename COT>
class AOBasisBuilder : public MPIObjectBase
{
public:
  enum
  {
    DONOT_EXPAND    = 0,
    GAUSSIAN_EXPAND = 1,
    NATURAL_EXPAND,
    CARTESIAN_EXPAND,
    MOD_NATURAL_EXPAND,
    DIRAC_CARTESIAN_EXPAND
  };

private:
  bool addsignforM;
  int expandlm;
  std::string Morder;
  std::string sph;
  std::string basisType;
  std::string elementType;
  std::string Normalized;

  ///map for the radial orbitals
  std::map<std::string, int> RnlID;

  ///map for (n,l,m,s) to its quantum number index
  std::map<std::string, int> nlms_id;

public:
  AOBasisBuilder(const std::string& eName, Communicate* comm);

  bool put(xmlNodePtr cur);
  bool putH5(hdf_archive& hin);

  SPOSet* createSPOSetFromXML(xmlNodePtr cur) { return 0; }

  std::unique_ptr<COT> createAOSet(xmlNodePtr cur);
  std::unique_ptr<COT> createAOSetH5(hdf_archive& hin);

  int expandYlm(COT* aos, std::vector<int>& all_nl, int expandlm = DONOT_EXPAND);
};
} // namespace qmcplusplus
#endif
