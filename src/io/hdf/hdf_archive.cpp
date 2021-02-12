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


#include <config.h>
#include "hdf_archive.h"
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

void hdf_archive::set_access_plist(bool request_pio, Communicate* comm)
{
  access_id = H5P_DEFAULT;
  myComm    = comm;
  if (comm && comm->size() > 1) //for parallel communicator
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
      H5Pset_fapl_mpio(access_id, comm->getMPI(), info);
      xfer_plist = H5Pcreate(H5P_DATASET_XFER);
      // enable parallel collective I/O
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
      use_phdf5 = true;
#else
      use_phdf5 = false;
#endif
    }
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
void hdf_archive::set_access_plist(bool request_pio, boost::mpi3::communicator& comm)
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
      H5Pset_fapl_mpio(access_id, &comm, info);
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
  hid_t g = H5Gopen(p, aname.c_str());
  if (g < 0)
    return false;
  H5Gclose(g);
  return true;
}

hid_t hdf_archive::push(const std::string& gname, bool createit)
{
  if (Mode[NOIO] || file_id == is_closed)
    return is_closed;
  hid_t p = group_id.empty() ? file_id : group_id.top();
  hid_t g = H5Gopen(p, gname.c_str());
  if (g < 0 && createit)
  {
    g = H5Gcreate(p, gname.c_str(), 0);
  }
  if (g != is_closed)
    group_id.push(g);
  return g;
}

} // namespace qmcplusplus
