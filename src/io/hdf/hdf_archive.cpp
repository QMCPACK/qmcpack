//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//		      Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "hdf_archive.h"
#ifdef HAVE_MPI
#include "mpi3/communicator.hpp"
#endif
#include "Message/Communicate.h"

namespace qmcplusplus
{
hdf_archive::~hdf_archive()
{
#if defined(ENABLE_PHDF5)
  if (xfer_plist != H5P_DEFAULT)
    H5Pclose(xfer_plist);
  if (access_id != H5P_DEFAULT)
    H5Pclose(access_id);
#endif
  close();
  H5Eset_auto2(H5E_DEFAULT, err_func, client_data);
}

void hdf_archive::close()
{
  while (!group_id.empty())
  {
    hid_t gid = group_id.top();
    group_id.pop();
    H5Gclose(gid);
  }
  if (file_id != is_closed)
    H5Fclose(file_id);
  file_id = is_closed;
}

void hdf_archive::set_access_plist()
{
  access_id = H5P_DEFAULT;
  Mode.set(IS_PARALLEL, false);
  Mode.set(IS_MASTER, true);
  Mode.set(NOIO, false);
}

void hdf_archive::set_access_plist(Communicate* comm, bool request_pio)
{
  access_id = H5P_DEFAULT;
  if (comm && comm->size() > 1) //for parallel communicator
  {
    bool use_phdf5 = false;
#if defined(ENABLE_PHDF5)
    if (request_pio)
    {
      // enable parallel I/O
      MPI_Info info = MPI_INFO_NULL;
      access_id     = H5Pcreate(H5P_FILE_ACCESS);
#if H5_VERSION_GE(1, 10, 0)
      H5Pset_all_coll_metadata_ops(access_id, true);
      H5Pset_coll_metadata_write(access_id, true);
#endif
      H5Pset_fapl_mpio(access_id, comm->getMPI(), info);
      xfer_plist = H5Pcreate(H5P_DATASET_XFER);
      // enable parallel collective I/O
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
      use_phdf5 = true;
    }
#endif
    Mode.set(IS_PARALLEL, use_phdf5);
    Mode.set(IS_MASTER, !comm->rank());
    if (request_pio && !use_phdf5)
      Mode.set(NOIO, comm->rank()); // master only
    else
      Mode.set(NOIO, false); // pio or all.
  }
  else
  {
    Mode.set(IS_PARALLEL, false);
    Mode.set(IS_MASTER, true);
    Mode.set(NOIO, false);
  }
}

#ifdef HAVE_MPI
void hdf_archive::set_access_plist(boost::mpi3::communicator& comm, bool request_pio)
{
  access_id = H5P_DEFAULT;
  if (comm.size() > 1) //for parallel communicator
  {
    bool use_phdf5 = false;
    if (request_pio)
    {
#if defined(ENABLE_PHDF5)
      // enable parallel I/O
      MPI_Info info = MPI_INFO_NULL;
      access_id     = H5Pcreate(H5P_FILE_ACCESS);
#if H5_VERSION_GE(1, 10, 0)
      H5Pset_all_coll_metadata_ops(access_id, true);
      H5Pset_coll_metadata_write(access_id, true);
#endif
      H5Pset_fapl_mpio(access_id, comm.get(), info);
      xfer_plist = H5Pcreate(H5P_DATASET_XFER);
      // enable parallel collective I/O
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
      use_phdf5 = true;
#else
      use_phdf5 = false;
#endif
    }
    Mode.set(IS_PARALLEL, use_phdf5);
    Mode.set(IS_MASTER, !comm.rank());
    if (request_pio && !use_phdf5)
      Mode.set(NOIO, comm.rank()); // master only
    else
      Mode.set(NOIO, false); // pio or all.
  }
  else
  {
    Mode.set(IS_PARALLEL, false);
    Mode.set(IS_MASTER, true);
    Mode.set(NOIO, false);
  }
}
#endif


bool hdf_archive::create(const std::string& fname, unsigned flags)
{
  if (Mode[NOIO])
    return true;
  if (!(Mode[IS_PARALLEL] || Mode[IS_MASTER]))
    throw std::runtime_error("Only create file in parallel or by master but not every rank!");
  close();
  file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access_id);
  return file_id != is_closed;
}

bool hdf_archive::open(const std::string& fname, unsigned flags)
{
  if (Mode[NOIO])
    return true;
  close();
  file_id = H5Fopen(fname.c_str(), flags, access_id);
  return file_id != is_closed;
}

bool hdf_archive::is_group(const std::string& aname)
{
  if (Mode[NOIO])
    return true;
  if (file_id == is_closed)
    return false;
  hid_t p = group_id.empty() ? file_id : group_id.top();
  p       = (aname[0] == '/') ? file_id : p;

  if (H5Lexists(p, aname.c_str(), H5P_DEFAULT) > 0)
  {
    H5O_info_t oinfo;
    oinfo.type = H5O_TYPE_UNKNOWN;

#if H5_VERSION_GE(1, 12, 0)
    H5Oget_info_by_name3(p, aname.c_str(), &oinfo, H5O_INFO_BASIC, H5P_DEFAULT);
#else
    H5Oget_info_by_name(p, aname.c_str(), &oinfo, H5P_DEFAULT);
#endif

    if (oinfo.type != H5O_TYPE_GROUP)
      return false;
    return true;
  }
  else
  {
    return false;
  }
}

hid_t hdf_archive::push(const std::string& gname, bool createit)
{
  hid_t g = is_closed;
  if (Mode[NOIO] || file_id == is_closed)
    return is_closed;
  hid_t p = group_id.empty() ? file_id : group_id.top();

  H5O_info_t oinfo;
  oinfo.type = H5O_TYPE_UNKNOWN;
  if (H5Lexists(p, gname.c_str(), H5P_DEFAULT) > 0)
  {
#if H5_VERSION_GE(1, 12, 0)
    H5Oget_info_by_name3(p, gname.c_str(), &oinfo, H5O_INFO_BASIC, H5P_DEFAULT);
#else
    H5Oget_info_by_name(p, gname.c_str(), &oinfo, H5P_DEFAULT);
#endif
  }

  if ((oinfo.type != H5O_TYPE_GROUP) && createit)
  {
    g = H5Gcreate2(p, gname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  else
  {
    g = H5Gopen2(p, gname.c_str(), H5P_DEFAULT);
  }
  if (g != is_closed)
    group_id.push(g);
  return g;
}

} // namespace qmcplusplus
