//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//		      Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "hdf_archive.h"
#include <H5Ipublic.h>
#include <stdexcept>
#ifdef HAVE_MPI
#include "mpi3/communicator.hpp"
#endif
#include "Message/Communicate.h"

namespace qmcplusplus
{

hdf_error_suppression hide_hdf_errors;

hdf_archive::~hdf_archive()
{
  close();
  H5Pclose(xfer_plist);
  H5Pclose(file_apl_);
  H5Pclose(file_cpl_);
  H5Pclose(lcpl_id);
}

void hdf_archive::close()
{
  while (!group_id.empty())
  {
    hid_t gid = group_id.top();
    group_id.pop();
    group_names.pop_back();
    H5Gclose(gid);
  }
  if (file_id != is_closed)
    H5Fclose(file_id);
  file_id = is_closed;
}

void hdf_archive::create_basic_plist()
{
  file_apl_ = H5Pcreate(H5P_FILE_ACCESS);
  if (file_apl_ == H5I_INVALID_HID)
    throw std::runtime_error("hdf_archive failed to create file access properties!");
  file_cpl_ = H5Pcreate(H5P_FILE_CREATE);
  if (file_cpl_ == H5I_INVALID_HID)
  {
    H5Pclose(file_apl_);
    throw std::runtime_error("hdf_archive failed to create file create properties!");
  }
  lcpl_id = H5Pcreate(H5P_LINK_CREATE);
  if (lcpl_id == H5I_INVALID_HID)
  {
    H5Pclose(file_apl_);
    H5Pclose(file_cpl_);
    throw std::runtime_error("hdf_archive failed to create link create properties!");
  }
  if (H5Pset_create_intermediate_group(lcpl_id, true) < 0)
  {
    H5Pclose(file_apl_);
    H5Pclose(file_cpl_);
    H5Pclose(lcpl_id);
    throw std::runtime_error("hdf_archive failed to set link create properties!");
  }
  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  if (xfer_plist == H5I_INVALID_HID)
  {
    H5Pclose(file_apl_);
    H5Pclose(file_cpl_);
    H5Pclose(lcpl_id);
    throw std::runtime_error("hdf_archive failed to create dataset xfer properties!");
  }
}

void hdf_archive::set_access_plist(Communicate& comm, bool request_pio)
{
#ifdef HAVE_MPI
  set_access_plist(comm.getCommMPI3(), request_pio);
#else
  if (request_pio)
    throw std::runtime_error("hdf_archive collective parallel I/O requested but ENABLE_PHDF5 is not enabled!");
  create_basic_plist();
#endif
}

#ifdef HAVE_MPI
void hdf_archive::set_access_plist(boost::mpi3::communicator& comm, bool request_pio)
{
#if !defined(ENABLE_PHDF5)
  if (request_pio)
    throw std::runtime_error("hdf_archive collective parallel I/O requested but ENABLE_PHDF5 is not enabled!");
#endif
  create_basic_plist();
#if defined(ENABLE_PHDF5)
  if (comm.size() > 1) //for parallel communicator
  {
    // enable parallel I/O
    MPI_Info info = MPI_INFO_NULL;
    if (H5Pset_fapl_mpio(file_apl_, comm.get(), info) < 0 ||
        (request_pio &&
         (H5Pset_all_coll_metadata_ops(file_apl_, true) < 0 || H5Pset_coll_metadata_write(file_apl_, true) < 0 ||
          H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)))
    {
      H5Pclose(file_apl_);
      H5Pclose(file_cpl_);
      H5Pclose(lcpl_id);
      H5Pclose(xfer_plist);
      throw std::runtime_error("hdf_archive failed to set parallel HDF5 properties!");
    }
  }
#endif
}
#endif

bool hdf_archive::create(const std::filesystem::path& fname, unsigned mode_flags)
{
  close();
  possible_filename_ = fname;
  if (read_only_predetermined_)
    throw std::runtime_error("Cannot create file in read-only mode!");
  file_id = H5Fcreate(fname.c_str(), mode_flags, file_cpl_, file_apl_);
  return file_id != is_closed;
}

bool hdf_archive::open(const std::filesystem::path& fname, unsigned mode_flags)
{
  close();
  possible_filename_ = fname;
  if (read_only_predetermined_ && mode_flags != H5F_ACC_RDONLY)
    throw std::runtime_error("Only H5F_ACC_RDONLY is allowed in read-only mode!");
  file_id = H5Fopen(fname.c_str(), mode_flags, file_apl_);
  return file_id != is_closed;
}

bool hdf_archive::is_group(const std::string& aname)
{
  if (file_id == is_closed)
    return false;
  hid_t p = group_id.empty() ? file_id : group_id.top();
  p       = (aname[0] == '/') ? file_id : p;

  if (H5Lexists(p, aname.c_str(), H5P_DEFAULT) > 0)
  {
#if H5_VERSION_GE(1, 12, 0)
    H5O_info2_t oinfo;
#else
    H5O_info_t oinfo;
#endif
    oinfo.type = H5O_TYPE_UNKNOWN;
    herr_t status;
#if H5_VERSION_GE(1, 12, 0)
    status = H5Oget_info_by_name3(p, aname.c_str(), &oinfo, H5O_INFO_BASIC, H5P_DEFAULT);
#else
    status = H5Oget_info_by_name(p, aname.c_str(), &oinfo, H5P_DEFAULT);
#endif

    if (status < 0 || oinfo.type != H5O_TYPE_GROUP)
      return false;
    return true;
  }
  else
  {
    return false;
  }
}

void hdf_archive::push(const std::string& gname, bool createit)
{
  hid_t g = is_closed;

  if (file_id == is_closed)
    throw std::runtime_error("Failed to open group \"" + gname +
                             "\" because file is not open.  Expected file: " + possible_filename_);

  hid_t p = group_id.empty() ? file_id : group_id.top();

#if H5_VERSION_GE(1, 12, 0)
  H5O_info2_t oinfo;
#else
  H5O_info_t oinfo;
#endif
  oinfo.type = H5O_TYPE_UNKNOWN;
  if (H5Lexists(p, gname.c_str(), H5P_DEFAULT) > 0)
  {
    herr_t status;
#if H5_VERSION_GE(1, 12, 0)
    status = H5Oget_info_by_name3(p, gname.c_str(), &oinfo, H5O_INFO_BASIC, H5P_DEFAULT);
#else
    status = H5Oget_info_by_name(p, gname.c_str(), &oinfo, H5P_DEFAULT);
#endif
    if (status < 0)
      oinfo.type = H5O_TYPE_UNKNOWN;
  }

  if ((oinfo.type != H5O_TYPE_GROUP) && createit)
  {
    if (read_only_predetermined_)
      throw std::runtime_error("Cannot create group in read-only mode!");
    g = H5Gcreate2(p, gname.c_str(), lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
  }
  else
  {
    g = H5Gopen2(p, gname.c_str(), H5P_DEFAULT);
  }
  if (g != is_closed)
  {
    group_id.push(g);
    group_names.push_back(gname);
  }

  if (!createit && g < 0)
    throw std::runtime_error("Group \"" + gname + "\" not found in file " + possible_filename_ +
                             ". Group path: " + group_path_as_string());
}

void hdf_archive::push(const hdf_path& gname, bool createit) { push(gname.string(), createit); }

std::string hdf_archive::group_path_as_string() const
{
  std::string group_path;
  for (auto it = group_names.begin(); it != group_names.end(); it++)
  {
    group_path += *it;
    if (it != group_names.end() - 1)
      group_path += "/";
  }

  return group_path;
}

} // namespace qmcplusplus
