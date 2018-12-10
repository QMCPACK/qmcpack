//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


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
#if defined(HAVE_LIBBOOST)
#if !defined(__bgq__)
#include <io/hdf_boost_smvector.h>
#endif
#endif
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
  enum {IS_PARALLEL=0, IS_MASTER, NOIO};
  static const hid_t is_closed=-1;
  /** bitset of the io mode
   * Mode[IS_PARALLEL] : true, if parallel
   * Mode[IS_MASTER] : true, if the node is master
   * Mode[NOIO] : true, if I/O is not performed
   */
  std::bitset<4> Mode;
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
  ///Pointer to communicator
  Communicate* myComm;
  ///FILO to handle H5Group
  std::stack<hid_t> group_id;
  /** constructor
   * @param c communicator
   * @param request_pio turns on parallel I/O,
   *        if true and PHDF5 is available, hdf_archive is in parallel collective IO mode
   *        if true and PHDF5 is not available, hdf_archive is in master-only IO mode
   *        if false, hdf_archive is in independent IO mode
   */
  hdf_archive(Communicate* c=nullptr, bool request_pio=false);
  ///destructor
  ~hdf_archive();

  ///set the access property
  void set_access_plist(bool request_pio, Communicate* comm);

  ///return true if parallel i/o
  inline bool is_parallel() const
  {
    return Mode[IS_PARALLEL];
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
  bool is_group(const std::string& aname);

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
    if(Mode[NOIO]) return true;
    if(!(Mode[IS_PARALLEL]||Mode[IS_MASTER])) std::runtime_error("Only write data in parallel or by master but not every rank!");
    hid_t p=group_id.empty()? file_id:group_id.top();
    h5data_proxy<T> e(data);
    return e.write(p,aname,xfer_plist);
  }

  template<typename T> bool read(T& data, const std::string& aname)
  {
    if(Mode[NOIO]) return true;
    hid_t p=group_id.empty()? file_id:group_id.top();
    h5data_proxy<T> e(data);
    return e.read(p,aname,xfer_plist);
  }

  inline void unlink(const std::string& aname)
  {
    if(Mode[NOIO]) return;
    hid_t p=group_id.empty()? file_id:group_id.top();
    herr_t status=H5Gunlink(p,aname.c_str());
  }
};
}
#endif
