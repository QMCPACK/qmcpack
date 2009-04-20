//////////////////////////////////////////////////////////////////
// (c) Copyright 2009- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
/** @file profiler.h
 * @brief interface to profiler libraries
 */
#ifndef QMCPLUSPLUS_PROFILER_H
#define QMCPLUSPLUS_PROFILER_H

#if defined(PROFILING_ON)
#include <TAU.h>
#else
//empty macros
#define TAU_PROFILE(a,b,c) 
#define TAU_INIT(argc,argv) 
#endif
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3727 $   $Date: 2009-04-03 11:49:36 -0400 (Fri, 03 Apr 2009) $
 * $Id: profiler.h 3727 2009-04-03 15:49:36Z jnkim $ 
 ***************************************************************************/
