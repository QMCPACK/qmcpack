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
  RecordArray(size_t nentries, size_t nparams) : Base(nentries, nparams) {}

  /// Change the size
  void resize(size_t nentries, size_t nparams)
  {
    Base::resize(nentries, nparams);
  }

  int getNumOfEntries() const { return Base::rows(); }
  int getNumOfParams() const { return Base::cols(); }
};
} // namespace qmcplusplus

#endif
