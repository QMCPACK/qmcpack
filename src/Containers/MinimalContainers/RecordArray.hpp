//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_RECORDARRAY_H
#define QMCPLUSPLUS_RECORDARRAY_H

#include <stdio.h>
#include <vector>

// This is a data structure with two indices.  The data is a list of samples,
// and each sample has an entry for multiple parameters.
// Think of a spreadsheet where the samples are the rows and the parameters are the columns.
// This could be implemented with nested vectors (std::vector<std::vector<T>> records),
// but that requires specifying the order of access everywhere it is used.  Any
// changes to the storage or indexing order require changes everywhere.  This data
// structure is an indirection between the specification of the access and the
// underlying storage.
//
// In principle the underlying data could be switched from AoS (record-oriented)
// to SoA (column-oriented), but the performance of accessing this structure is
// unlikely to matter enough to warrant such functionality.
//
// The calls to setValue and getValue use the parameter index first, and the
// sample entry index second.  It would nice to create different types for these
// indices so the correct ordering could be enforced by the compiler.
//

namespace qmcplusplus
{
template<typename T>
class RecordArray
{
public:
  void setValue(int param_idx, int entry_idx, T value)
  {
    storage_[param_idx * param_dim_ + entry_idx * nparam_] = value;
  }

  T getValue(int param_idx, int entry_idx) const
  {
    return storage_[param_idx * param_dim_ + entry_idx * entry_dim_];
  }

  RecordArray() : nparam_(0), nentry_(0) {}

  /// Constructor specifying the number of parameters and entries

  RecordArray(int nparam, int nentry) : nparam_(nparam), nentry_(nentry)
  {
    storage_.resize(nparam_ * nentry_);
    // AoS (or record-oriented) - param is fastest varying
    param_dim_ = 1;
    entry_dim_ = nparam_;
  }

  /// Change the size
  void resize(int nparam, int nentry)
  {
    nparam_ = nparam;
    nentry_ = nentry;

    param_dim_ = 1;
    entry_dim_ = nparam_;

    storage_.resize(nparam_ * nentry_);
  }

  int nparam() const { return nparam_;}

  int nentry() const { return nentry_;}

private:
  // Number of parameters (columns)
  int nparam_;
  // Number of entries (rows)
  int nentry_;

  // These control the ordering of the storage
  int param_dim_;
  int entry_dim_;

  std::vector<T> storage_;
};
} // namespace qmcplusplus

#endif
