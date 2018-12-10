//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file profiler.h
 * @brief interface to profiler libraries and event id
 *
 * Experimental. Do not include in the header file.
 */
#ifndef QMCPLUSPLUS_PROFILER_H
#define QMCPLUSPLUS_PROFILER_H

#include <config.h>
#include <stack>
#if defined(VTRACE)
#include <vt_user.h>
namespace qmcplusplus {
  struct Profile
  {
    std::stack<std::string> domains;
    inline Profile(const char* name, const char* tag,  unsigned id)
    {
      domains.push(name);
      VT_USER_START(name);
    }
    inline ~Profile()
    {
      pop();
    }
    inline void push(const char* name)
    {
      domains.push(name);
      VT_USER_START(name);
    }

    inline void pop()
    {
      std::string name=domains.top();
      VT_USER_END(name.c_str());
      domains.pop();
    }
  };
}
#elif defined(PROFILEING_ON)
#include <TAU.h>
namespace qmcplusplus {
  struct Profile
  {
    unsigned ID;
    inline Profile(const char* name, const char* tag,  unsigned id)
      :ID(id)
    {
      TAU_PROFILE(name,tag,id);
    }
    inline ~Profile()
    { }
    inline void push(const char* name)
    { }
    inline void pop()
    { }
  };
}
//#elif defined(PROFILEING_HPCT_ON)
//#include <libhpc.h>
#else
namespace qmcplusplus {
  struct Profile
  {
    inline Profile(const char* name, const char* tag,  unsigned id)
    { }
    inline ~Profile()
    { }
    inline void push(const char* name)
    { }

    inline void pop()
    { }
  };
}
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

#define QMC_WFS_0_EVENT 0x0300
#define QMC_WFS_1_EVENT 0x0310
#define QMC_WFS_2_EVENT 0x0320
#define QMC_WFS_3_EVENT 0x0330
#define QMC_WFS_4_EVENT 0x0340
#define QMC_WFS_5_EVENT 0x0350
#define QMC_WFS_6_EVENT 0x0360
#define QMC_WFS_7_EVENT 0x0370
#define QMC_WFS_8_EVENT 0x0380

#endif
