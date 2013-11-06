//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#include <QMCWaveFunctions/SPOInfo.h>
#include <limits>


namespace qmcplusplus
{

  typedef QMCTraits::RealType RealType;

  const int      SPOInfo::no_index       = -1;
  const int      SPOInfo::no_degeneracy  = -1;
  const RealType SPOInfo::no_energy      = 1e99;

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

  void SPOInfo::report(const string& pad)
  {
    if(has_index())
      app_log()<<pad<<"index      = "<< index      <<endl;
    else
      app_log()<<pad<<"index      = not assigned"  <<endl;
    if(has_energy())
      app_log()<<pad<<"energy     = "<< energy     <<endl;
    else
      app_log()<<pad<<"energy     = not assigned"  <<endl;
    if(has_degeneracy())
      app_log()<<pad<<"degeneracy = "<< degeneracy <<endl;
    else
      app_log()<<pad<<"degeneracy = not assigned"  <<endl;
    app_log().flush();
  }
}
