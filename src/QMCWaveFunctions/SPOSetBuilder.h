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

#include <memory>
#include <vector>
#include <string>
#include "Message/MPIObjectBase.h"
#include "QMCWaveFunctions/SPOSetInfo.h"
#include "QMCWaveFunctions/SPOSetInputInfo.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "hdf/hdf_archive.h"

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
class SPOSetBuilder : public QMCTraits, public MPIObjectBase
{
public:
  using SPOPool_t  = std::map<std::string, SPOSet*>;
  using indices_t  = std::vector<int>;
  using energies_t = std::vector<RealType>;

  /// whether implementation conforms only to legacy standard
  bool legacy;

  /// state info of all possible states available in the basis
  std::vector<std::unique_ptr<SPOSetInfo>> states;

  SPOSetBuilder(const std::string& type_name, Communicate* comm);
  virtual ~SPOSetBuilder() {}

  /// reserve space for states (usually only one set, multiple for e.g. spin dependent einspline)
  void reserve_states(int nsets = 1);

  /// allow modification of state information
  inline void modify_states(int index = 0) { states[index]->modify(); }

  /// clear state information
  inline void clear_states(int index = 0) { states[index]->clear(); }

  /// create an sposet from xml and save the resulting SPOSet
  std::unique_ptr<SPOSet> createSPOSet(xmlNodePtr cur);

  const std::string& getTypeName() const { return type_name_; }

protected:
  /// create an sposet from xml (legacy)
  virtual std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) = 0;

  /// create an sposet from a general xml request
  virtual std::unique_ptr<SPOSet> createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input_info);

  /// type name of the SPO objects built by this builder.
  const std::string type_name_;
};

} // namespace qmcplusplus
#endif
