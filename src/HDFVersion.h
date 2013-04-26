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

#include <config.h>
#include <io/hdf.h>

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

struct HDFVersion//: public HDFAttribIOBase
{
  //enumeration to get version value
  enum {MAJOR=0, MINOR};
  typedef TinyVector<int,2> data_type;
  data_type version;

  inline HDFVersion():
    version(QMCPLUSPLUS_VERSION_MAJOR,QMCPLUSPLUS_VERSION_MINOR)
  { }

  inline explicit HDFVersion(int m, int n):version(m,n)
  { }

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

  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    h5data_proxy<data_type> vin(version);
    return vin.read(grp,aname,xfer_plist);
  }

  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    h5data_proxy<data_type> vout(version);
    return vout.write(grp,aname,xfer_plist);
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

/** specialize h5data_proxy for HDFVersion */
template<>
struct h5data_proxy<HDFVersion>
{
  HDFVersion& ref;
  h5data_proxy(HDFVersion& a):ref(a) {}
  inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    h5data_proxy<HDFVersion::data_type> vin(ref.version);
    return vin.read(grp,aname,xfer_plist);
  }
  inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
  {
    h5data_proxy<HDFVersion::data_type> vout(ref.version);
    return vout.write(grp,aname,xfer_plist);
  }
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1575 $   $Date: 2007-01-04 09:36:37 -0600 (Thu, 04 Jan 2007) $
 * $Id: HDFWalkerIO.h 1575 2007-01-04 15:36:37Z jnkim $
 ***************************************************************************/
