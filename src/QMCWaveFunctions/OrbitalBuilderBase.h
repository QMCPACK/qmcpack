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
    
    
/**@file OrbitalBuilderBase.h
 *@brief declaration of the base class for many-body wavefunction.
 */
#ifndef QMCPLUSPLUS_TRIALORBITALBUILDERBASE_H
#define QMCPLUSPLUS_TRIALORBITALBUILDERBASE_H

#include "QMCWaveFunctions/TrialWaveFunction.h"
#include <map>

/**@defgroup WFSBuilder Orbital builder group
 * @brief Builder classes to add OrbitalComponent to a TrialWaveFunction
 */
namespace qmcplusplus
{

/**@ingroup WFSBuilder
 * @brief An abstract class for wave function builders
 */
class OrbitalBuilderBase: public MPIObjectBase
{

public:

  typedef TrialWaveFunction::RealType  RealType;
  typedef TrialWaveFunction::ValueType ValueType;
  typedef TrialWaveFunction::PosType   PosType;
  typedef TrialWaveFunction::GradType  GradType;
  typedef std::map<std::string,ParticleSet*> PtclPoolType;

  /////level of printing
  //static int print_level;

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
  /// the element name for a finite-difference linear response wavefunction
  static std::string fdlrwfn_tag;
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
  /// the element name for the basis set: basis set contains multiple basis elements
  static std::string basisset_tag;
  /// the element name for a group of basis functions: basis contains multiple basisfunc elements
  static std::string basis_tag;
  /// the element name for a basis function
  static std::string basisfunc_tag;
  /// the element name for an ion wavefunction
  static std::string ionorb_tag;
  /// the element name for a backflow transformation
  static std::string backflow_tag;
  /// the element name for a multi slater determinant wavefunction
  static std::string multisd_tag;
  /// the element name for a SPO scanner
  static std::string sposcanner_tag;

  //@}

  /** constructor
   * @param p target particle set
   * @param psi target many-body wavefunction
   *
   * Each builder class adds an object to compose a many-body wavefunction
   * targetPsi. The position of targetPtcl is related to targetPsi's
   * capability to return a value and derivatives \f$\Psi[\{R\}]\f$ .
   */
  OrbitalBuilderBase(ParticleSet& p, TrialWaveFunction& psi);

  virtual ~OrbitalBuilderBase();
  /// process a xml node at cur
  virtual bool put(xmlNodePtr cur) = 0;

protected:
  /// reference to the particle set on which targetPsi is defined
  ParticleSet& targetPtcl;

  /// reference to the many-body wavefucntion to which each derived class add a term
  TrialWaveFunction& targetPsi;

  /// xmlNode operated by this object
  xmlNodePtr myNode;

  /// child builder
  std::vector<OrbitalBuilderBase*> Children;
};

}
#endif
