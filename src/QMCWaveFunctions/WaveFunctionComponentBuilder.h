//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file WaveFunctionComponentBuilder.h
 *@brief declaration of the base class for many-body wavefunction.
 */
#ifndef QMCPLUSPLUS_TRIALORBITALBUILDERBASE_H
#define QMCPLUSPLUS_TRIALORBITALBUILDERBASE_H

#include <map>
#include <memory>
#include "Message/MPIObjectBase.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"

/**@defgroup WFSBuilder Orbital builder group
 * @brief Builder classes to build WaveFunctionComponent
 */
namespace qmcplusplus
{
/**@ingroup WFCBuilder
 * @brief An abstract class for wave function builders
 */
class WaveFunctionComponentBuilder : public MPIObjectBase
{
public:
  using RealType  = WaveFunctionComponent::RealType;
  using ValueType = WaveFunctionComponent::ValueType;
  using PosType   = WaveFunctionComponent::PosType;
  using GradType  = WaveFunctionComponent::GradType;
  using PSetMap   = std::map<std::string, std::unique_ptr<ParticleSet>>;

  /** \ingroup XMLTags
   *@{
   *@brief reserved tags for the elements associated with the many-body wavefunctions
   */
  /// the element name for a wavefunction
  static std::string wfs_tag;
  /// the element name for a parameter
  static std::string param_tag;
  /// the element name for a distancetable
  static std::string dtable_tag;
  /// the element name for jatrow
  static std::string jastrow_tag;
  /// the element name for a set of Slater determinants, contains 1..* Slater determinants
  static std::string detset_tag;
  /// the element name for a Slater determinant, contains 1..* determinants
  static std::string sd_tag;
  /// the element name for a determinant, may contain (0..*) orbital or parameter element
  static std::string det_tag;
  /// the element name for a released node determinant, may contain (0..*) orbital or parameter element
  static std::string rn_tag;
  /// the element name for single-particle orbital
  static std::string spo_tag;
  /// the element name for single-particle orbital set
  static std::string sposet_tag;
  /// the element name for an ion wavefunction
  static std::string ionorb_tag;
  /// the element name for a backflow transformation
  static std::string backflow_tag;
  /// the element name for a multi slater determinant wavefunction
  static std::string multisd_tag;

  //@}

  /** constructor
   * @param comm communicator
   * @param p target particle set
   *
   * Each builder class builds an object for composing a many-body wavefunction.
   */
  WaveFunctionComponentBuilder(Communicate* comm, ParticleSet& p) : MPIObjectBase(comm), targetPtcl(p), myNode(NULL) {}

  virtual ~WaveFunctionComponentBuilder() = default;
  /// process a xml node at cur
  virtual std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) = 0;

protected:
  /// reference to the particle set on which targetPsi is defined
  ParticleSet& targetPtcl;

  /// xmlNode operated by this object
  xmlNodePtr myNode;
};

} // namespace qmcplusplus
#endif
