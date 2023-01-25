//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Steven Hahn, hahnse@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Steven Hahn, hahnse@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF5_PATH_H
#define QMCPLUSPLUS_HDF5_PATH_H

#include <initializer_list>
#include <string>
#include <string_view>

namespace qmcplusplus
{
class hdf_path
{
public:
  /** Constructs an empty path.
   */
  hdf_path();
  /** Constructs a path whose pathname is the same as that of p
   *  @param p pathname  
   */
  hdf_path(std::string_view p);
  /** Appends elements to the path with a directory separator
   *  @param p pathname
   */
  hdf_path& append(std::string_view p);
  hdf_path& append(const hdf_path& p);
  /** Concatenates two paths without introducing a directory separator
   *  @param p pathname
   *  @return *this
   */
  hdf_path& concat(std::string_view p);
  hdf_path& concat(const hdf_path& p);
  /** Appends elements to the path with a directory separator
   *  @param p pathname
   *  @return *this
   */
  hdf_path& operator/=(std::string_view p);
  hdf_path& operator/=(const hdf_path& p);
  /** Concatenates two paths without introducing a directory separator
   *  @param p pathname
   *  @return *this
   */
  hdf_path& operator+=(std::string_view p);
  hdf_path& operator+=(const hdf_path& p);
  /** Remove last subgroup in path
   *  @return *this
   */
  hdf_path& remove_subgroup();
  /** Replace the last subgroup in path
    * @param p subgroup to replace last in path
    * @return *this
    */
  hdf_path& replace_subgroup(std::string_view p);
  /** Return std::string for use with HDF5
    * @return path as a std::string
    */
  const std::string& string() const { return path_; }
  /** Does the path absolute (true) or relative to the current group in hdf_archive (false)?
   *  @return does the path begin with '/' 
   */
  bool has_root_directory() const;

private:
  std::string path_;
};

/** concatenates two paths with a directory separator
 * @param lhs left-hand side of directory separator
 * @param rhs right-hand side of directory separator
 * @return path
 */
hdf_path operator/(const hdf_path& lhs, const hdf_path& rhs);
hdf_path operator/(const hdf_path& lhs, const std::string& rhs);
hdf_path operator/(const hdf_path& lhs, const char* rhs);
/** Checks whether lhs and rhs are equal.
 * @param lhs first path
 * @param rhs second path
 * @return result
 */
bool operator==(const hdf_path& lhs, const hdf_path& rhs) noexcept;

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_HDF5_PATH_H
