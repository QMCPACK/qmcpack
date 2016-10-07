//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include <QMCWaveFunctions/SPOInfo.h>
#include <limits>


namespace qmcplusplus
{

  typedef QMCTraits::RealType RealType;

  const int      SPOInfo::no_index       = -1;
  const int      SPOInfo::no_degeneracy  = -1;
  const RealType SPOInfo::no_energy      = std::numeric_limits<RealType>::max();

  SPOInfo::SPOInfo()
  {
    index      = no_index; 
    degeneracy = no_degeneracy;
    energy     = no_energy;
  }

  SPOInfo::SPOInfo(int orb_index,RealType en)
  {
    index      = orb_index;
    degeneracy = no_degeneracy;
    energy     = en;
  }

  void SPOInfo::report(const std::string& pad)
  {
    if(has_index())
      app_log()<<pad<<"index      = "<< index      << std::endl;
    else
      app_log()<<pad<<"index      = not assigned"  << std::endl;
    if(has_energy())
      app_log()<<pad<<"energy     = "<< energy     << std::endl;
    else
      app_log()<<pad<<"energy     = not assigned"  << std::endl;
    if(has_degeneracy())
      app_log()<<pad<<"degeneracy = "<< degeneracy << std::endl;
    else
      app_log()<<pad<<"degeneracy = not assigned"  << std::endl;
    app_log().flush();
  }
}
