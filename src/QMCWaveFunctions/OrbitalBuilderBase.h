//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  typedef map<string,ParticleSet*> PtclPoolType;

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
  vector<OrbitalBuilderBase*> Children;
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
