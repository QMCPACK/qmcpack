////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by:
// Fionn Malone, malone14@llnl.gov, Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef HDF5_CONSISTENCY_HELPER_HPP
#define HDF5_CONSISTENCY_HELPER_HPP

#include <vector>

#include "hdf/hdf_multi.h"
#include "hdf/hdf_archive.h"

#include "AFQMC/config.0.h"


namespace qmcplusplus
{
namespace afqmc
{
// Helper functions for reading integral data from HDF5 files.
// Currently only read one-dimensional integrals into std::vector<ValueType> and two-dimensional
// integrals into boost::multi::array<ValueType,2>, so only handle these overloads.

/** read data from filespace (name) to buffer (out).
 * If QMC_COMPLEX = 1 and T = RealType then handle this through copy.
 * If QMC_COMPLEX = 0 and T = ComplexType then abort with error message.
 * @return true if successful. false indiciates name not found in dump.
 */
template<typename T>
bool readComplexOrReal(hdf_archive& dump, std::string name, std::vector<T>& out)
{
  std::vector<int> shape;
  int ndim = 1; // vector
  if (!dump.getShape<std::vector<T>>(name, shape))
  {
    // name not found in dump.
    return false;
  }
#ifdef QMC_COMPLEX
  if (shape.size() == ndim + 1)
  {
    // Nothing.
    dump.readEntry(out, name);
    return true;
  }
  else if (shape.size() == ndim)
  {
    // Real integrals.
    std::vector<RealType> out_real(out.size());
    dump.readEntry(out_real, name);
    std::copy_n(out_real.begin(), out_real.size(), out.begin());
    return true;
  }
  else
  {
    app_log() << " Error reading " << name << " dataspace. Shape mismatch.\n";
    APP_ABORT("");
    return false;
  }
#else
  if (shape.size() == ndim + 1)
  {
    app_log() << " Error: Found complex integrals with QMC_COMPLEX=0.\n";
    APP_ABORT(" Please recompile with QMC_COMPLEX=1 or generate real integrals if appropriate.\n");
    return false;
  }
  else
  {
    dump.readEntry(out, name);
    return true;
  }
#endif
}

template<typename T>
bool readComplexOrReal(hdf_archive& dump, std::string name, boost::multi::array<T, 2>& out)
{
  std::vector<int> shape;
  int ndim = 2; // matrix
  if (!dump.getShape<boost::multi::array<T, 2>>(name, shape))
  {
    // name not found in dump.
    return false;
  }
#ifdef QMC_COMPLEX
  if (shape.size() == ndim + 1)
  {
    // Nothing.
    dump.readEntry(out, name);
    return true;
  }
  else if (shape.size() == ndim)
  {
    // Real integrals.
    boost::multi::array<RealType, 2> out_real({shape[0], shape[1]});
    dump.readEntry(out_real, name);
    std::copy_n(out_real.origin(), out_real.num_elements(), out.origin());
    return true;
  }
  else
  {
    app_log() << " Error reading " << name << " dataspace. Shape mismatch.\n";
    APP_ABORT("");
    return false;
  }
#else
  if (shape.size() == ndim + 1)
  {
    app_log() << " Error: Found complex integrals with QMC_COMPLEX=0.\n";
    APP_ABORT(" Please recompile with QMC_COMPLEX=1 or generate real integrals if appropriate.\n");
    return false;
  }
  else
  {
    dump.readEntry(out, name);
    return true;
  }
#endif
}

} // namespace afqmc
} // namespace qmcplusplus

#endif
