//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_RECORDARRAY_H
#define QMCPLUSPLUS_RECORDARRAY_H

#include "OhmmsPETE/OhmmsMatrix.h"

// This is a data structure with two indices. The data is a list of samples,
// and each sample has an entry for multiple parameters.
// Think of a spreadsheet where the samples are the rows and the parameters are the columns.
//
namespace qmcplusplus
{
template<typename T>
class RecordArray : public Matrix<T>
{
public:
  using Base = Matrix<T>;

  RecordArray() {}

  /// Constructor specifying the number of parameters and entries

  RecordArray(int nparam, int nentry) : Base(nentry, nparam), num_entries_(nentry), num_params_(nparam) {}

  /// Change the size
  void resize(int nparam, int nentry)
  {
    num_entries_ = nentry;
    num_params_  = nparam;
    Base::resize(nentry, nparam);
  }

  int getNumOfParams() const { return num_params_; }

  int getNumOfEntries() const { return num_entries_; }

private:
  // Number of entries (rows)
  int num_entries_ = 0;
  // Number of parameters (columns)
  int num_params_ = 0;
};
} // namespace qmcplusplus

#endif
