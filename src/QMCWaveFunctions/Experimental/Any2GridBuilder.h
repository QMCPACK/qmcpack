//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_ANY2RADIALGRIDFUNCTOR_H
#define QMCPLUSPLUS_ANY2RADIALGRIDFUNCTOR_H

#include "QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.h"

namespace qmcplusplus
{

/**Class to convert SlaterTypeOrbital to a radial orbital on a log grid.
 *
 * For a center,
 *   - only one grid is used
 *   - any number of radial orbitals
 */
struct Any2GridBuilder: public RGFBuilderBase
{

  ///constructor
  Any2GridBuilder(xmlNodePtr cur=NULL);

  ///implement the virtual function
  bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);
  bool putCommon(xmlNodePtr cur);

  bool Normalized;
  RealType m_rcut;
  QuantumNumberType m_nlms;

  void addGaussian(xmlNodePtr cur);
  void addSlater(xmlNodePtr cur);
  void addPade(xmlNodePtr cur);
};

}
#endif
