//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef OHMMS_ATOMICHARTREEFOCK_TYPES_H
#define OHMMS_ATOMICHARTREEFOCK_TYPES_H

#include "AtomicHF/YlmRnlSet.h"

/**@namespace ohmmshf
 *@brief Define basic data types for the applications.
 * In order to reduce complier-time complexity and to enable switching
 * between C++ libraries for array and expression template, 
 * basic data types are defined.
 */
namespace ohmmshf {
  typedef YlmRnlSet<double> HFAtomicOrbitals;
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
