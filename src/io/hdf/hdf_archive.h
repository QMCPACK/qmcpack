//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF5_ARCHIVE_H
#define QMCPLUSPLUS_HDF5_ARCHIVE_H

#include <config.h>
#include "hdf_datatype.h"
#include "hdf_dataspace.h"
#include "hdf_dataproxy.h"
#include "hdf_error_suppression.h"
#include "hdf_path.h"
#include "hdf_pete.h"
#include "hdf_stl.h"
#include "hdf_hyperslab.h"

#include <bitset>
#include <filesystem>
#include <stack>

#ifdef HAVE_MPI
namespace boost
{
namespace mpi3
{
class communicator;
}
} // namespace boost
#endif

class Communicate;

namespace qmcplusplus
{
/** class to handle hdf file
 */
class hdf_archive
{
private:
  enum
  {
    IS_PARALLEL = 0,
    IS_MASTER,
    NOIO
  };
  static const hid_t is_closed = -1;
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
  /// Link creation property list identifier
  hid_t lcpl_id;
  ///FILO to handle H5Group
  std::stack<hid_t> group_id;
  ///Track group names corresponding to group_id
  std::vector<std::string> group_names;

  /** Name of file that hdf_archive thinks is open.
   *  This may not correspond to the actual file because the open call failed,
   *  or the file was closed. This information is useful for debugging.
   */
  std::string possible_filename_;

  ///set the access property
  void set_access_plist(Communicate* comm, bool request_pio);
#ifdef HAVE_MPI
  void set_access_plist(boost::mpi3::communicator& comm, bool request_pio);
#endif
  void set_access_plist();

public:
  /** constructor
   * @param c communicator
   * @param request_pio turns on parallel I/O,
   *        if true and PHDF5 is available, hdf_archive is in parallel collective IO mode
   *        if true and PHDF5 is not available, hdf_archive is in master-only IO mode
   *        if false, hdf_archive is in independent IO mode
   */
  template<class Comm = Communicate*>
  hdf_archive(Comm c, bool request_pio = false)
      : file_id(is_closed), access_id(H5P_DEFAULT), xfer_plist(H5P_DEFAULT), lcpl_id(H5P_DEFAULT)
  {
    if (!hdf_error_suppression::enabled)
      throw std::runtime_error("HDF5 library warnings and errors not suppressed from output.\n");
    set_access_plist(c, request_pio);
  }

  hdf_archive() : file_id(is_closed), access_id(H5P_DEFAULT), xfer_plist(H5P_DEFAULT), lcpl_id(H5P_DEFAULT)
  {
    if (!hdf_error_suppression::enabled)
      throw std::runtime_error("HDF5 library warnings and errors not suppressed from output.\n");
    set_access_plist();
  }
  ///destructor
  ~hdf_archive();

  ///return true if parallel i/o
  inline bool is_parallel() const { return Mode[IS_PARALLEL]; }

  ///return true if master in parallel i/o
  inline bool is_master() const { return Mode[IS_MASTER]; }

  ///return file_id. should be only be used for connecting to old codes when porting
  hid_t getFileID() const { return file_id; }

  /** create a file
   * @param fname name of hdf5 file
   * @param flags i/o mode
   * @return true, if creation is successful
   */
  bool create(const std::filesystem::path& fname, unsigned flags = H5F_ACC_TRUNC);

  /** open a file
   * @param fname name of hdf5 file
   * @param flags i/o mode
   * @return file_id, if open is successful
   */
  bool open(const std::filesystem::path& fname, unsigned flags = H5F_ACC_RDWR);

  ///close all the open groups and file
  void close();

  ///flush a file
  inline void flush()
  {
    if (file_id != is_closed)
      H5Fflush(file_id, H5F_SCOPE_LOCAL);
  }

  ///return true if the file is closed
  inline bool closed() { return file_id == is_closed; }

  /** check if aname is a group
   * @param aname group's name
   * @return true, if aname exists and it is a group
   */
  bool is_group(const std::string& aname);

  /** check if aname is a dataset
   * @param aname dataset's name
   * @return true, if aname exists and it is a dataset
   */
  bool is_dataset(const std::string& aname)
  {
    if (Mode[NOIO])
      return true;
    hid_t p = group_id.empty() ? file_id : group_id.top();
    int dummy_data;
    h5data_proxy<int> e(dummy_data);
    return e.check_existence(p, aname);
  }

  /** check if aname is a dataset of type T
   * @param aname group's name
   * @return true, if aname is a dataset of type T
   */
  template<typename T>
  bool is_dataset_of_type(const std::string& aname)
  {
    if (Mode[NOIO])
      return true;
    hid_t p = group_id.empty() ? file_id : group_id.top();
    T dummy_data;
    h5data_proxy<T> e(dummy_data);
    return e.check_type(p, aname);
  }

  /** return the top of the group stack
   */
  inline hid_t top() const { return group_id.empty() ? file_id : group_id.top(); }

  /** check if any groups are open
   *  group stack will have entries if so
   *  @return true if any groups are open
   */
  inline bool open_groups() { return group_id.empty(); }

  /** push a group to the group stack
   * @param gname name of the group
   * @param createit if true, group is create when missing
   */
  void push(const std::string& gname, bool createit = true);
  void push(const hdf_path& gname, bool createit = true);


  inline void pop()
  {
    if (file_id == is_closed || group_id.empty())
      return;
    hid_t g = group_id.top();
    group_id.pop();
    group_names.pop_back();
    herr_t err = H5Gclose(g);
    if (err < 0)
      throw std::runtime_error("H5Gclose failed with error.");
  }

  /** Return a string representation of the current group stack
   */
  std::string group_path_as_string() const;

  /** read the shape of multidimensional filespace from the group aname
   * this function can be used to query dataset for preparing containers.
   * The dimensions contributed by T is excluded.
   * See how exactly user dimensions are calculated in getDataShape function definition.
   * @return true if successful
   */
  template<typename T>
  bool getShape(const std::string& aname, std::vector<int>& sizes_out)
  {
    if (Mode[NOIO])
      return true;
    hid_t p = group_id.empty() ? file_id : group_id.top();
    return getDataShape<T>(p, aname, sizes_out);
  }

  /** write the data to the group aname and return status
   * use write() for inbuilt error checking
   * @return true if successful
   */
  template<typename T>
  bool writeEntry(T& data, const std::string& aname)
  {
    if (Mode[NOIO])
      return true;
    if (!(Mode[IS_PARALLEL] || Mode[IS_MASTER]))
      throw std::runtime_error("Only write data in parallel or by master but not every rank!");
    hid_t p = group_id.empty() ? file_id : group_id.top();
    h5data_proxy<typename std::remove_const<T>::type> e(data);
    return e.write(data, p, aname, xfer_plist);
  }

  /** write the data to the group aname and check status
   * runtime error is issued on I/O error
   */
  template<typename T>
  void write(T& data, const std::string& aname)
  {
    if (!writeEntry(data, aname))
    {
      throw std::runtime_error("HDF5 write failure in hdf_archive::write " + aname);
    }
  }

  /** write the container data with a specific shape and check status
   * @param data container, linear storage required.
   * @param shape shape on the hdf file
   * @param aname dataset name in the file
   * runtime error is issued on I/O error
   */
  template<typename T, typename IT, std::size_t RANK>
  void writeSlabReshaped(T& data, const std::array<IT, RANK>& shape, const std::string& aname)
  {
    std::array<hsize_t, RANK> globals, counts, offsets;
    for (int dim = 0; dim < RANK; dim++)
    {
      globals[dim] = static_cast<hsize_t>(shape[dim]);
      counts[dim]  = static_cast<hsize_t>(shape[dim]);
      offsets[dim] = 0;
    }

    hyperslab_proxy<T, RANK> pxy(data, globals, counts, offsets);
    write(pxy, aname);
  }

  /** read the data from the group aname and return status
   * use read() for inbuilt error checking
   * @return true if successful
   */
  template<typename T, typename = std::enable_if_t<!std::is_const<T>::value>>
  bool readEntry(T& data, const std::string& aname)
  {
    if (Mode[NOIO])
      return true;
    hid_t p = group_id.empty() ? file_id : group_id.top();
    h5data_proxy<T> e(data);
    return e.read(data, p, aname, xfer_plist);
  }

  /** read the data from the group aname and check status
   * runtime error is issued on I/O error
   */
  template<typename T, typename = std::enable_if_t<!std::is_const<T>::value>>
  void read(T& data, const std::string& aname)
  {
    if (!readEntry(data, aname))
    {
      throw std::runtime_error("HDF5 read failure in hdf_archive::read " + aname);
    }
  }

  /** read file dataset with a specific shape into a container and check status
   * @param data container, linear storage required.
   * @param shape shape on the hdf file
   * @param aname dataset name in the file
   * runtime error is issued on I/O error
   */
  template<typename T, typename IT, std::size_t RANK, typename = std::enable_if_t<!std::is_const<T>::value>>
  void readSlabReshaped(T& data, const std::array<IT, RANK>& shape, const std::string& aname)
  {
    std::array<hsize_t, RANK> globals, counts, offsets;
    for (int dim = 0; dim < RANK; dim++)
    {
      globals[dim] = static_cast<hsize_t>(shape[dim]);
      counts[dim]  = static_cast<hsize_t>(shape[dim]);
      offsets[dim] = 0;
    }

    hyperslab_proxy<T, RANK> pxy(data, globals, counts, offsets);
    read(pxy, aname);
  }

  /** read a portion of the data from the group aname and check status
   * runtime error is issued on I/O error
   *
   * note the readSpec array must have dimensionality corresponding to the dataset,
   * values for a dimension must be [0,num_entries-1] for that dimension to specify
   * which value to hold and a -1 to grab all elements from that dimension
   * for example, if the dataset was [5,2,6] and the vector contained (2,1,-1),
   * this would grab 6 elements corresponding to [2,1,:]
   */
  template<typename T, typename IT, std::size_t RANK, typename = std::enable_if_t<!std::is_const<T>::value>>
  void readSlabSelection(T& data, const std::array<IT, RANK>& readSpec, const std::string& aname)
  {
    std::array<hsize_t, RANK> globals, counts, offsets;
    for (int dim = 0; dim < RANK; dim++)
    {
      globals[dim] = 0;
      if (readSpec[dim] < 0)
      {
        counts[dim]  = 0;
        offsets[dim] = 0;
      }
      else
      {
        counts[dim]  = 1;
        offsets[dim] = static_cast<hsize_t>(readSpec[dim]);
      }
    }

    hyperslab_proxy<T, RANK> pxy(data, globals, counts, offsets);
    read(pxy, aname);
  }

  inline void unlink(const std::string& aname)
  {
    if (Mode[NOIO])
      return;
    hid_t p       = group_id.empty() ? file_id : group_id.top();
    herr_t status = H5Ldelete(p, aname.c_str(), H5P_DEFAULT);
  }
};

} // namespace qmcplusplus
#endif
