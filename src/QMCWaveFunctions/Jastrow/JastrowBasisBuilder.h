//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
  JastrowBasisBuilder(ParticleSet& els, ParticleSet& ions, const std::string& functype, bool usespline=false);
  ///destructor
  ~JastrowBasisBuilder();

  bool put(xmlNodePtr cur);

  /** do not create anything. */
  SPOSetBase* createSPOSetFromXML(xmlNodePtr cur)
  {
    return 0;
  }
  ///size of blocks
  std::vector<int> SizeOfBasisPerCenter;
private:
  ///target ParticleSet
  ParticleSet& targetPtcl;
  ///source ParticleSet
  ParticleSet& sourcePtcl;
  ///boolean to choose numerical or analytic form
  bool UseSpline;
  ///function name
  std::string FuncType;
  ///save AtomiBasisBuilder<RFB>*
  std::map<std::string,BasisSetBuilder*> aoBuilders;
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
  void printRadialFunctors(const std::string& elementType, COT* aoBasis);
};
}
#endif
