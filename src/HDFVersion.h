//////////////////////////////////////////////////////////////////
// (c) Copyright 2007- by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_HDFWALKERVERSION_H
#define QMCPLUSPLUS_HDFWALKERVERSION_H

#include "config.h"
#include "Numerics/HDFNumericAttrib.h"
#define ENABLE_PHDFLIB
namespace qmcplusplus
{

  namespace hdf
  {
    //1st level
    const char version[]="version";
    const char main_state[]="state_0";
    const char config_group[]="config_collection";

    //2nd level for main_state
    const char random[]="random_state";
    const char walkers[]="walkers";
    const char num_walkers[]="number_of_walkers";
    const char energy_history[]="energy_history";
    const char norm_history[]="norm_history";

    //2nd level for config_group
    const char num_blocks[]="NumOfConfigurations";
    const char append_walkers[]="config_";

    //unused
    const char coord[]="coord";
  }

  struct HDFVersion: public HDFAttribIOBase
  {
    //enumeration to get version value
    enum {MAJOR=0, MINOR};
    typedef TinyVector<int,2> data_type;
    data_type version;

    inline HDFVersion(): version(QMCPLUSPLUS_VERSION_MAJOR,QMCPLUSPLUS_VERSION_MINOR){}
    inline explicit HDFVersion(int m, int n):version(m,n){}

    inline void write(hid_t gid, const char* name)
    {
      HDFAttribIO<data_type> h(version);
      h.write(gid,name);
    }

    inline void read(hid_t gid, const char* name)
    {
      HDFAttribIO<data_type> h(version);
      h.read(gid,name);
    }

    inline int operator[](int i) const
    {
      return version[i];
    }

    //could be general to D 
    inline bool operator==(const HDFVersion &other) const 
    {
      return (version[0] == other.version[0] && version[1] == other.version[1]);
    }

    inline bool operator!=(const HDFVersion &other) const 
    {
      return (version[0] != other.version[0] || version[1] != other.version[1]);
    }

    ///limited to 100 for each version field
    inline int serialized() const 
    {
      return version[0]*100+version[1];
    }

    inline bool operator>=(const HDFVersion &other) const
    {
      return serialized()>=other.serialized();
    }
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1575 $   $Date: 2007-01-04 09:36:37 -0600 (Thu, 04 Jan 2007) $
 * $Id: HDFWalkerIO.h 1575 2007-01-04 15:36:37Z jnkim $ 
 ***************************************************************************/
