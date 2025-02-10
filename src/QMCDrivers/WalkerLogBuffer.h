//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_WALKERLOGBUFFER_H
#define QMCPLUSPLUS_WALKERLOGBUFFER_H


#include <Configuration.h>
#include <unordered_set>
#include "OhmmsPETE/OhmmsArray.h"
#include "hdf/hdf_archive.h"
#include "Utilities/TimerManager.h"

namespace qmcplusplus
{

/// Basic data types for walker log data
struct WLog
{
  using Int  = long;                 // integer type
  using Real = OHMMS_PRECISION_FULL; // real type
  using Comp = std::complex<Real>;   // complex type
#ifndef QMC_COMPLEX
  using PsiVal = Real; // wavefunction type
#else
  using PsiVal = Comp; // wavefunction type
#endif
};


/** Record for an individual walker quantity being logd.
 *
 *    Helper struct for WalkerLogBuffer.
 */
struct WalkerQuantityInfo
{
  /// quantity name
  std::string name;
  /// support up to 4D array quantity
  enum
  {
    D0 = 0,
    D1,
    D2,
    D3,
    DMAX
  };
  /// array dimension
  size_t dimension;
  /// total size
  size_t size;
  /// size of 1 unit of data
  size_t unit_size;
  /// array shape
  TinyVector<size_t, DMAX> shape;
  /// starting row index in buffer
  size_t buffer_start;
  /// end range in buffer row
  size_t buffer_end;

  WalkerQuantityInfo(const std::string& name_,
                     size_t unit_size_,
                     size_t buffer_start_,
                     size_t n1 = 1,
                     size_t n2 = 0,
                     size_t n3 = 0,
                     size_t n4 = 0)
  {
    name         = name_;
    unit_size    = unit_size_;
    buffer_start = buffer_start_;
    shape[D0]    = n1;
    shape[D1]    = n2;
    shape[D2]    = n3;
    shape[D3]    = n4;

    dimension = 0;
    size      = 1;
    for (size_t d = 0; d < DMAX; ++d)
      if (shape[d] > 0)
      {
        dimension++;
        size *= shape[d];
      }
    buffer_end = buffer_start + size * unit_size;
  }
};


/** Data buffer for walker log quantities.
 *
 *    Each row in the buffer contains all quantities for one walker from a single step.
 *    Rows are added throughout an MC block.
 *    See WalkerLogCollector::collect()
 *
 *    Buffer data is written to HDF at the end of each MC block.
 *    See WalkerLogManager::writeBuffers()
 */
template<typename T>
class WalkerLogBuffer
{
public:
  /// label for this data in HDF file
  std::string label;
  /// marks first WalkerLogCollector::collect() call
  bool first_collect;
  /// HDF file pointer
  hsize_t hdf_file_pointer;

private:
  static constexpr std::string_view my_name_{"WalkerLogBuffer"};
  enum Timer
  {
    COLLECT = 0,
    ADD_ROW,
    WRITE,
    RESET,
    MAKE_NEW_ROW
  };
  static constexpr std::array<std::string_view, 5> suffixes_{"collect", "add_row", "write", "reset", "make_new_row"};
  static TimerNameList_t<Timer> create_names(const std::string_view& my_name)
  {
    TimerNameList_t<Timer> timer_names;
    using namespace std::string_literals;
    std::string prefix{"WalkerLog:"s + std::string{my_name} + "::"s};
    for (std::size_t i = 0; i < suffixes_.size(); ++i)
      timer_names.push_back({static_cast<Timer>(i), prefix + std::string{suffixes_[i]}});
    return timer_names;
  }
  TimerList_t walker_log_buffer_timers_;

  /// index of current quantity during WalkerLogCollector::collect()
  size_t quantity_index;
  /** buffer row location data for each walker quantity
   *    used to populate "data_layout" field in HDF file
   */
  std::vector<WalkerQuantityInfo> quantity_info;
  /// total size of walker data stored in a buffer row
  size_t walker_data_size;
  /// the walker data buffer itself
  Array<T, 2> buffer;

public:
  WalkerLogBuffer();

  /// current number of rows in the data buffer
  inline size_t nrows() { return buffer.size(0); }

  /// current number of columns in the data buffer (row size)
  inline size_t ncols() { return buffer.size(1); }

  /// resize the buffer to zero
  inline void resetBuffer() { buffer.resize(0, buffer.size(1)); }

  /// reset member variables at end of each WalkerLogCollector::collect() call
  inline void resetCollect()
  {
    if (quantity_index != quantity_info.size())
      throw std::runtime_error(
          "WalkerLogBuffer quantity_index has not been moved through all quantities prior during collect() call.");
    first_collect  = false;
    quantity_index = 0;
  }

  /// compare row size of this buffer to another one
  inline bool sameAs(const WalkerLogBuffer<T>& ref) { return buffer.size(1) == ref.buffer.size(1); }

  /// collect data for a single walker quantity of scalar type into the current buffer row
  void collect(const std::string& name, const T& value);

  /** collect data for a single walker quantity of array type into the current buffer row
   *  Only D = 1,2 are actually used by any dependent code so only those D's are explicit instantiated
   */
  template<unsigned D>
  void collect(const std::string& name, Array<T, D> arr);

  /** collect data for a single walker quantity of complex array type into the current buffer row
   *  Only D = 1,2 are actually used by any dependent code so only those D's are explicit instantiated
   */
  template<unsigned D>
  void collect(const std::string& name, Array<std::complex<T>, D> arr);

  /// add a data row from another buffer to this one
  void addRow(WalkerLogBuffer<T> other, size_t i);

  /// write a summary of quantities in the buffer
  void writeSummary(std::string pad = "  ");

  /// write the data_layout for all walker quantities in the HDF file
  void registerHDFData(hdf_archive& f);

  /// write the buffer data into the HDF file
  void writeHDF(hdf_archive& f);

  /// write the buffer data into the HDF file
  void writeHDF(hdf_archive& f, hsize_t& file_pointer);

private:
  /// make space as quantities are added to the buffer for the first time
  void resetRowSize(size_t row_size);

  /// allocate a full new row at the end of the buffer
  void makeNewRow();
};

// explicit instantiations
extern template class WalkerLogBuffer<WLog::Int>;
extern template class WalkerLogBuffer<WLog::Real>;
// and for templated member functions.
extern template void WalkerLogBuffer<WLog::Int>::collect<2>(const std::string& name, Array<WLog::Int, 2>);
extern template void WalkerLogBuffer<WLog::Real>::collect<2>(const std::string& name, Array<WLog::Real, 2>);
extern template void WalkerLogBuffer<WLog::Real>::collect<1>(const std::string& name, Array<WLog::Real, 1>);
extern template void WalkerLogBuffer<WLog::Real>::collect<2>(const std::string& name,
                                                             Array<std::complex<WLog::Real>, 2>);
extern template void WalkerLogBuffer<WLog::Real>::collect<1>(const std::string& name,
                                                             Array<std::complex<WLog::Real>, 1>);

} // namespace qmcplusplus

#endif
