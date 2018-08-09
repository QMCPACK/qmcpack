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


///\file BsplineDeviceCUDA.h
#ifndef QMCPLUSPLUS_BSPLINEDEVICECUDA_H
#define QMCPLUSPLUS_BSPLINEDEVICECUDA_H

#include <iostream>
#include "QMCWaveFunctions/BsplineFactory/BsplineDevice.h"
#include "einspline/bspline_structs.h"
namespace qmcplusplus
{

/** 
 * \class BsplineDeviceCUDA
 */
template<typename ST, unsigned D>
class BsplineDeviceCUDA : BsplineDevice<BsplineDeviceCUDA<ST,D>, ST, D>
{
public:
  //So on the host side we still use
  using SingleBsplineType=UBspline_3d_d;
  
  void implementation()
  {
    std::cout<< "implemented for CUDA\n";
  }
};

}
#endif
