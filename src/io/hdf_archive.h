//////////////////////////////////////////////////////////////////
// (c) Copyright 2010- by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_HDF5_ARCHIVE_H
#define QMCPLUSPLUS_HDF5_ARCHIVE_H

#include <Configuration.h>
#include <io/hdf_datatype.h>
#include <io/hdf_dataspace.h>
#include <io/hdf_dataproxy.h>
#if defined(HAVE_LIBHDF5)
#include <io/hdf_pete.h>
#include <io/hdf_stl.h>
#include <io/hdf_hyperslab.h>
#endif
#include <stack>
#include <bitset>

class Communicate;

namespace qmcplusplus
{
/** class to handle hdf file
 */
struct hdf_archive
{
  enum {IS_PARALLEL=0, NOIO};
  static const hid_t is_closed=-1;
  /** bitset of the io mode
   * Mode[IS_PARALLEL] : true, if collective
   * Mode[NOIO] : true, if I/O is not performed
   */
  bitset<4> Mode;
  ///file id
  hid_t file_id;
  ///access id
  hid_t access_id;
  ///transfer property
  hid_t xfer_plist;
  ///error type
  H5E_auto_t err_func;
  ///error handling
  void *client_data;
  ///FILO to handle H5Group
  std::stack<hid_t> group_id;
  /** constructor
   * @param c communicator
   * @param use_collective turn on/off collective
   */
  hdf_archive(Communicate* c=0, bool use_collective=false);
  ///destructor
  ~hdf_archive();

  ///set the access property
  void set_access_plist(bool use_collective, Communicate* comm);

  ///return true if collective i/o
  inline bool is_collective() const
  {
    return Mode[0];
  }

  /** create a file
   * @param fname name of hdf5 file
   * @param flags i/o mode
   * @return true, if creation is successful
   */
  bool create(const std::string& fname, unsigned flags=H5F_ACC_TRUNC);

  /** open a file
   * @param fname name of hdf5 file
   * @param flags i/o mode
   * @return file_id, if open is successful
   */
  bool open(const std::string& fname,unsigned flags=H5F_ACC_RDWR);

  ///close all the open groups and file
  void close();

  ///flush a file
  inline void flush()
  {
    if(file_id!=is_closed)
      H5Fflush(file_id,H5F_SCOPE_LOCAL);
  }

  /** check if aname is a group
   * @param aname group's name
   * @return true, if aname exists and it is a group
   */
  bool is_group(const string& aname);

  /** return the top of the group stack
   */
  inline hid_t top() const
  {
    return group_id.empty()?is_closed:group_id.top();
  }

  /** push a group to the group stack
   * @param gname name of the group
   * @param createit if true, group is create when missing
   */
  hid_t push(const std::string& gname, bool createit=true);

  inline void pop()
  {
    if(file_id==is_closed || group_id.empty())
      return;
    hid_t g=group_id.top();
    group_id.pop();
    H5Gclose(g);
  }

  template<typename T> bool write(T& data, const std::string& aname)
  {
    if(Mode[NOIO])
      return true;
    hid_t p=group_id.empty()? file_id:group_id.top();
    h5data_proxy<T> e(data);
    return e.write(p,aname,xfer_plist);
  }

  template<typename T> bool read(T& data, const std::string& aname)
  {
    if(Mode[NOIO])
      return true;
    hid_t p=group_id.empty()? file_id:group_id.top();
    h5data_proxy<T> e(data);
    return e.read(p,aname,xfer_plist);
  }

  inline void unlink(const std::string& aname)
  {
    if(Mode[NOIO])
      return;
    hid_t p=group_id.empty()? file_id:group_id.top();
    herr_t status=H5Gunlink(p,aname.c_str());
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 894 $   $Date: 2006-02-03 10:52:38 -0600 (Fri, 03 Feb 2006) $
 * $Id: hdf.h 894 2006-02-03 16:52:38Z jnkim $
 ***************************************************************************/
