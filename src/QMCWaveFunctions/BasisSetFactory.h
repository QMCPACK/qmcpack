//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_BASISSETFACTORY_H
#define QMCPLUSPLUS_BASISSETFACTORY_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/BasisSetBase.h"

namespace qmcplusplus
{

  ///writes info about contained sposets to stdout
  void write_basis_builders(const string& pad="");

  /**returns a named sposet from the global pool
   *  only use in serial portion of execution
   *  ie during initialization prior to threaded code
   */
  SPOSetBase* get_sposet(const string& name);



/** derived class from OrbitalBuilderBase
 */
class BasisSetFactory: public OrbitalBuilderBase
{

public:

  ///set of basis set builders resolved by type
  static map<string,BasisSetBuilder*> basis_builders;

  /** constructor
   * \param els reference to the electrons
   * \param psi reference to the wavefunction
   * \param ions reference to the ions
   */
  BasisSetFactory(ParticleSet& els, TrialWaveFunction& psi, PtclPoolType& psets);

  ~BasisSetFactory();
  bool put(xmlNodePtr cur);

  BasisSetBuilder* createBasisSet(xmlNodePtr cur, xmlNodePtr rootNode=NULL);

  SPOSetBase* createSPOSet(xmlNodePtr cur);

  void build_sposet_collection(xmlNodePtr cur);

private:

  ///store the last builder, use if type not provided
  static BasisSetBuilder* last_builder;

  ///reference to the particle pool
  PtclPoolType& ptclPool;
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
