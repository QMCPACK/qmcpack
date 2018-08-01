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
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file BplineDevice.h
 *
 * Base class for BsplineDevice<DEVICE> 
 */

#ifndef QMCPLUSPLUS_BSPLINEDEVICE_H
#define QMCPLUSPLUS_BSPLINEDEVICE_H


namespace qmcplusplus
{

/** base class of any BsplineDevice
 *
 * Defines interface for storage of einspline
 */
template<class DEVICETYPE, typename ST, unsigned D>
class BsplineDevice
{
  void interface()
  {
    static_cast<DEVICETYPE*>(this)->implementation();
  }
};

}

#endif
