//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file SPOSetBuilder.h
 * @brief Declaration of a base class of SPOSet Builders
 */
#ifndef QMCPLUSPLUS_SPOSET_BUILDER_H
#define QMCPLUSPLUS_SPOSET_BUILDER_H

#include "Message/MPIObjectBase.h"
#include "QMCWaveFunctions/SPOSetInfo.h"
#include <QMCWaveFunctions/SPOSetInputInfo.h>
#include "QMCWaveFunctions/SPOSet.h"


namespace qmcplusplus
{
/** base class for the real SPOSet builder
 *
 * \warning {
 * We have not quite figured out how to use real/complex efficiently.
 * There are three cases we have to deal with
 * - real basis functions and real coefficients
 * - real basis functions and complex coefficients
 * - complex basis functions and complex coefficients
 * For now, we decide to keep both real and complex basis sets and expect
 * the user classes {\bf KNOW} what they need to use.
 * }
 */
struct SPOSetBuilder : public QMCTraits, public MPIObjectBase
{
  typedef std::map<std::string, SPOSet*> SPOPool_t;
  typedef std::vector<int> indices_t;
  typedef std::vector<RealType> energies_t;


  /// whether implementation conforms only to legacy standard
  bool legacy;

  /// state info of all possible states available in the basis
  std::vector<SPOSetInfo*> states;

  /// list of all sposets created by this builder
  std::vector<SPOSet*> sposets;

  SPOSetBuilder(Communicate* comm);
  virtual ~SPOSetBuilder() {}
  /// load from XML if there is a basisset
  virtual void loadBasisSetFromXML(xmlNodePtr cur) {}

  /// reserve space for states (usually only one set, multiple for e.g. spin dependent einspline)
  void reserve_states(int nsets = 1);

  /// allow modification of state information
  inline void modify_states(int index = 0) { states[index]->modify(); }

  /// clear state information
  inline void clear_states(int index = 0) { states[index]->clear(); }

  /// create an sposet from xml (legacy)
  virtual SPOSet* createSPOSetFromXML(xmlNodePtr cur) = 0;

  /// create an sposet from a general xml request
  virtual SPOSet* createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input_info);

  /// create an sposet from xml and save the resulting SPOSet
  SPOSet* createSPOSet(xmlNodePtr cur);
};

} // namespace qmcplusplus
#endif
