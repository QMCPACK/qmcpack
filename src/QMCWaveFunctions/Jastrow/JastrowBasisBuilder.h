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
#ifndef QMCPLUSPLUS_JASTROWBASISBUILDER_H
#define QMCPLUSPLUS_JASTROWBASISBUILDER_H

#include "QMCWaveFunctions/BasisSetBase.h"
//#include "QMCWaveFunctions/Jastrow/CBSOBuilder.h"
//#include "QMCWaveFunctions/LocalizedBasisSet.h"

namespace qmcplusplus
{


/** derived class from BasisSetBuilder
 *
 * Create a basis set of molecular orbital types as defined in MolecularOrbitalBasis
 * with radial wave functions on the radial grids.
 */
class JastrowBasisBuilder: public BasisSetBuilder
{
public:
  ///keep the current basis set
  BasisSetBase<RealType>* myBasisSet;
  /** constructor
   * \param els reference to the electrons
   * \param ions reference to the ions
   */
  JastrowBasisBuilder(ParticleSet& els, ParticleSet& ions, const string& functype, bool usespline=false);
  ///destructor
  ~JastrowBasisBuilder();

  bool put(xmlNodePtr cur);

  /** do not create anything. */
  SPOSetBase* createSPOSetFromXML(xmlNodePtr cur)
  {
    return 0;
  }
  ///size of blocks
  vector<int> SizeOfBasisPerCenter;
private:
  ///target ParticleSet
  ParticleSet& targetPtcl;
  ///source ParticleSet
  ParticleSet& sourcePtcl;
  ///boolean to choose numerical or analytic form
  bool UseSpline;
  ///function name
  string FuncType;
  ///save AtomiBasisBuilder<RFB>*
  map<string,BasisSetBuilder*> aoBuilders;
  /** create a localized basis set
   *
   * The template parameter RFBUILDER is a builder class for radial
   * functors.
   */
  template<typename RFBUILDER>
  void createLocalizedBasisSet(xmlNodePtr cur);

  /** print the basis set
   * @param elementType element name
   * @param aoBasis atomic orbtial basis associated with elementType
   */
  template<typename COT>
  void printRadialFunctors(const string& elementType, COT* aoBasis);
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1623 $   $Date: 2007-01-18 18:01:15 -0600 (Thu, 18 Jan 2007) $
 * $Id: JastrowBasisBuilder.h 1623 2007-01-19 00:01:15Z jnkim $
 ***************************************************************************/
