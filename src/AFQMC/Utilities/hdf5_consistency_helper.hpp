#ifndef HDF5_CONSISTENCY_HELPER_HPP
#define HDF5_CONSISTENCY_HELPER_HPP

#include <vector>

#include "io/hdf_multi.h"
#include "io/hdf_archive.h"

#include "AFQMC/config.0.h"


namespace qmcplusplus
{
namespace afqmc
{
// Helper functions for reading integral data from HDF5 files.
// Check if data matches type T, else copy if  QMC_COMPLEX=1, else abort.
template<typename T>
int readComplexOrReal(hdf_archive& dump, std::string name, std::vector<T>& out)
{
  std::vector<int> shape;
  int ndim = 1; // vector
  if(!dump.getShape<std::vector<T>>(name, shape)) {
    // name not found in dump.
    return false;
  }
#ifdef QMC_COMPLEX
  if(shape.size() == ndim+1) {
    // Nothing.
    dump.readEntry(out, name);
    return true;
  } else if (shape.size() == ndim) {
    // Real integrals.
    std::vector<RealType> out_real(out.size());
    dump.readEntry(out_real, name);
    std::copy_n(out_real.begin(), out_real.size(), out.begin());
    return true;
  } else {
    app_log() << " Error reading " << name << " dataspace. Shape mismatch.\n";
    APP_ABORT("");
  }
#else
  if(shape.size() == ndim+1) {
    app_log() << " Error: Found complex integrals with QMC_COMPLEX=0.\n";
    APP_ABORT(" Please recompile with QMC_COMPLEX=1 or generate real integrals if appropriate.\n");
  } else {
    dump.readEntry(out, name);
    return true;
  }
#endif
}

template<typename T>
bool readComplexOrReal(hdf_archive& dump, std::string name, boost::multi::array<T,2>& out)
{
  std::vector<int> shape;
  int ndim = 2; // matrix
  if(!dump.getShape<boost::multi::array<T,2>>(name, shape)) {
    // name not found in dump.
    return false;
  }
#ifdef QMC_COMPLEX
  if(shape.size() == ndim+1) {
    // Nothing.
    dump.readEntry(out, name);
    return true;
  } else if (shape.size() == ndim) {
    // Real integrals.
    boost::multi::array<RealType,2> out_real({shape[0],shape[1]});
    dump.readEntry(out_real, name);
    std::copy_n(out_real.origin(), out_real.num_elements(), out.origin());
    return true;
  } else {
    app_log() << " Error reading " << name << " dataspace. Shape mismatch.\n";
    APP_ABORT("");
  }
#else
  if(shape.size() == ndim+1) {
    app_log() << " Error: Found complex integrals with QMC_COMPLEX=0.\n";
    APP_ABORT(" Please recompile with QMC_COMPLEX=1 or generate real integrals if appropriate.\n");
  } else {
    dump.readEntry(out, name);
    return true;
  }
#endif
}

}
}

#endif
