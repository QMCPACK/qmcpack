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
 * @brief interface to profiler libraries and event id
 *
 * Experimental. Do not include in the header file.
 */
#ifndef QMCPLUSPLUS_PROFILER_H
#define QMCPLUSPLUS_PROFILER_H

#include <config.h>
//#define PROFILING_HPCT_ON

//provide minimal TAU interfaces
#if defined(PROFILING_ON)
#include <TAU.h>
#else
#define TAU_PROFILE(a,b,c)
#define TAU_INIT(argc,argv)
#endif

//provide minimal IBM HPCT interface
#if defined(PROFILING_HPCT_ON)
#include <libhpc.h>
#else
#define hpmInit(a,b)
#define hpmStart(a,b)
#define hpmStop(a)
#define hpmTerminate(a)
#endif

//QMC event IDS
#define QMC_MAIN_EVENT  0x0000

#define QMC_VMC_0_EVENT 0x0100
#define QMC_VMC_1_EVENT 0x0110
#define QMC_VMC_2_EVENT 0x0120
#define QMC_VMC_3_EVENT 0x0130
#define QMC_VMC_4_EVENT 0x0140
#define QMC_VMC_5_EVENT 0x0150
#define QMC_VMC_6_EVENT 0x0160
#define QMC_VMC_7_EVENT 0x0170
#define QMC_VMC_8_EVENT 0x0180

#define QMC_DMC_0_EVENT 0x0200
#define QMC_DMC_1_EVENT 0x0210
#define QMC_DMC_2_EVENT 0x0220
#define QMC_DMC_3_EVENT 0x0230
#define QMC_DMC_4_EVENT 0x0240
#define QMC_DMC_5_EVENT 0x0250
#define QMC_DMC_6_EVENT 0x0260
#define QMC_DMC_7_EVENT 0x0270
#define QMC_DMC_8_EVENT 0x0280

#define QMC_WFS_0_EVENT 0x1100
#define QMC_WFS_1_EVENT 0x1110
#define QMC_WFS_2_EVENT 0x1120
#define QMC_WFS_3_EVENT 0x1130
#define QMC_WFS_4_EVENT 0x1140
#define QMC_WFS_5_EVENT 0x1150
#define QMC_WFS_6_EVENT 0x1160
#define QMC_WFS_7_EVENT 0x1170
#define QMC_WFS_8_EVENT 0x1180

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3727 $   $Date: 2009-04-03 11:49:36 -0400 (Fri, 03 Apr 2009) $
 * $Id: profiler.h 3727 2009-04-03 15:49:36Z jnkim $
 ***************************************************************************/
