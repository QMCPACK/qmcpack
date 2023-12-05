//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Steven Hahn, hahnse@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Steven Hahn, hahnse@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "hdf_path.h"

namespace qmcplusplus
{
hdf_path::hdf_path() = default;

hdf_path::hdf_path(std::string_view p) : path_{p} {};

hdf_path& hdf_path::append(std::string_view p) { return append(hdf_path(p)); }

hdf_path& hdf_path::append(const hdf_path& p)
{
  if (p.has_root_directory())
    path_.clear();
  else if (!path_.empty() && path_.back() != '/')
    path_.push_back('/');
  return concat(p);
}

hdf_path& hdf_path::concat(std::string_view p)
{
  path_.append(p);
  return *this;
}

hdf_path& hdf_path::concat(const hdf_path& p) { return concat(p.string()); }

hdf_path& hdf_path::operator/=(std::string_view p) { return append(p); }

hdf_path& hdf_path::operator/=(const hdf_path& p) { return append(p); }

hdf_path& hdf_path::operator+=(std::string_view p) { return concat(p); }

hdf_path& hdf_path::operator+=(const hdf_path& p) { return concat(p); }

hdf_path& hdf_path::remove_subgroup()
{
  // if path == '/' do nothing
  if (path_.size() == 1 && has_root_directory())
    return *this;

  auto pos = path_.find_last_of('/');

  // don't eliminate root group
  if (pos == 0 && has_root_directory())
    pos++;

  // if path is relative and only one group is present, clear the string
  if (pos == std::string::npos)
    path_.clear();
  else
    path_.erase(pos, std::string::npos);
  return *this;
}

hdf_path& hdf_path::replace_subgroup(std::string_view p)
{
  remove_subgroup();
  return append(p);
}

bool hdf_path::has_root_directory() const { return (!path_.empty() && path_[0] == '/'); }

hdf_path operator/(const hdf_path& lhs, const hdf_path& rhs) { return hdf_path(lhs) /= rhs; }

hdf_path operator/(const hdf_path& lhs, const std::string& rhs) { return hdf_path(lhs) /= std::string_view(rhs); }

hdf_path operator/(const hdf_path& lhs, const char* rhs) { return hdf_path(lhs) /= std::string_view(rhs); }

bool operator==(const hdf_path& lhs, const hdf_path& rhs) noexcept { return lhs.string() == rhs.string(); }

} // namespace qmcplusplus
