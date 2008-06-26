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
namespace qmcplusplus
{

  namespace hdf
  {
    /// extension of a configuration file
    const char config_ext[]=".config.h5";

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
    const char qmc_status[]="qmc_status";

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
    hid_t h_xfer;

    inline HDFVersion():
    version(QMCPLUSPLUS_VERSION_MAJOR,QMCPLUSPLUS_VERSION_MINOR), 
    h_xfer(H5P_DEFAULT)
    { }

    inline explicit HDFVersion(int m, int n):version(m,n), h_xfer(H5P_DEFAULT)
    { }

    ~HDFVersion()
    {
      if(h_xfer != H5P_DEFAULT) H5Pclose(h_xfer);
    }

    inline void write(hid_t gid, const char* name)
    {
      HDFAttribIO<data_type> vin(version);
      vin.write(gid,name);
    }

    inline void read(hid_t gid, const char* name)
    {
      hid_t h1 = H5Dopen(gid,name);
      if(h1>-1)
      {
        hid_t ret = H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, h_xfer, version.begin());
        H5Dclose(h1);
      }
    }

    /** set parallel mode
     * @param parallel true, if the file is open by parallel
     */
    inline void setPmode(bool parallel)
    {
#if defined(H5_HAVE_PARALLEL)
      if(parallel)
      {
        h_xfer = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(h_xfer,H5FD_MPIO_INDEPENDENT);
      }
#endif
    }

    inline int operator[](int i) const
    {
      return version[i];
    }

    inline int& operator[](int i) 
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

    inline bool operator<(const HDFVersion &other) const
    {
      return serialized()<other.serialized();
    }
  };

  inline std::ostream& operator<<(std::ostream& out, const HDFVersion& v)
  {
    out << v[0] << " " << v[1];
    return out;
  }

  inline std::istream& operator>>(std::istream& is, HDFVersion& v)
  {
    is >> v[0] >> v[1];
    return is;
  }

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1575 $   $Date: 2007-01-04 09:36:37 -0600 (Thu, 04 Jan 2007) $
 * $Id: HDFWalkerIO.h 1575 2007-01-04 15:36:37Z jnkim $ 
 ***************************************************************************/
