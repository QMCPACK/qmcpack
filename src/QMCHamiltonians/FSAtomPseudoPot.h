//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and  Kenneth Esler
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
/** @file FSAtomPseudoPot.h
 * @brief Xml Parser Definition for FSAtom Pseudopotential Standard
 */
#ifndef QMCPLUSPLUS_FSATOMPSEDUOPOTENTIAL_H
#define QMCPLUSPLUS_FSATOMPSEDUOPOTENTIAL_H
#include "Configuration.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimLinearSpline.h"
#include "OhmmsData/AttributeSet.h"
namespace qmcplusplus
{

struct FSAtomPseudoPot: public QMCTraits
{
  typedef OneDimLinearSpline<RealType> return_type;
  typedef OneDimLinearSpline<RealType>::grid_type grid_type;
  return_type myFunc;
  int AngL;
  RealType Rcut;

  FSAtomPseudoPot(int l, RealType rc, grid_type* agrid):
    myFunc(agrid), AngL(l), Rcut(rc)
  { }

  ~FSAtomPseudoPot()
  { }

  void convert2RV()
  {
    for(int i=0; i<myFunc.size(); i++)
      myFunc(i) *= myFunc.r(i);
  }

  void convert2HartreeBohr(RealType sc, bool is_r_times_v);

  RealType getCutOff(RealType v0);

  /** create a LinearSpline<RealType>
   * @param sc scaling factor
   * @return a LinearSpline<RealType> on a LinearGrid
   */
  return_type* getLocalPot(RealType zeff);

  return_type* getNonLocalPot(FSAtomPseudoPot& vloc);

  bool put(xmlNodePtr cur);

};
} // namespace qmcPlusPlus
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
