//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_BSPLINEDEVICE_H
#define QMCPLUSPLUS_BSPLINEDEVICE_H

#include "spline2/MultiBspline.hpp"

namespace qmcplusplus
{

/** base class of any BsplineDevice
 *
 * Defines interface for storage of einspline
 */
template<class DEVICETYPE, typename ST, unsigned D>
class BsplineDevice
{
public:
  void interface()
  {
    static_cast<DEVICETYPE*>(this)->implementation();
  }

  
  void createSpline(MultiBspline<ST>& multi_bspline)
  {
    static_cast<DEVICETYPE*>(this)->createSpline_imp(multi_bspline);
  }

  void initDevice(MultiBspline<ST>& multi_bspline)
  {
    static_cast<DEVICETYPE*>(this)->initDevice_imp(multi_bspline);
  }

  void resizeStorage(size_t n, size_t nvals, int num_walkers)
  {
    static_cast<DEVICETYPE*>(this)->resizeStorage_imp(n, nvals, num_walkers);
  }
};

}

#endif
