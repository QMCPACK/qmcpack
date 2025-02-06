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

  /// collect data for a single walker quantity of array type into the current buffer row
  template<unsigned D>
  inline void collect(const std::string& name, Array<T, D> arr)
  {
    ScopedTimer timer(walker_log_buffer_timers_[Timer::COLLECT]);

    size_t n1 = arr.size(0);
    size_t n2, n3, n4;
    n2 = n3 = n4 = 0;
    if (D > 4)
      throw std::runtime_error("WalkerLogBuffer::collect  Only arrays up to dimension 4 are currently supported.");
    if (D > 1)
      n2 = arr.size(1);
    if (D > 2)
      n3 = arr.size(2);
    if (D > 3)
      n4 = arr.size(3);
    size_t irow = 0;
    if (first_collect)
    { // cache walker quantity info on first collect
      WalkerQuantityInfo wqi_(name, 1, walker_data_size, n1, n2, n3, n4);
      quantity_info.push_back(wqi_);
      walker_data_size = wqi_.buffer_end;
      resetRowSize(walker_data_size);
    }
    else
    { // make a new buffer row if needed
      if (quantity_index == 0)
        makeNewRow();
      irow = buffer.size(0) - 1;
    }
    // place the array walker quantity into the current buffer row
    auto& wqi   = quantity_info[quantity_index];
    auto& arr1d = arr.storage();
    for (size_t n = 0; n < arr1d.size(); ++n)
      buffer(irow, wqi.buffer_start + n) = arr1d[n];
    quantity_index++;
  }


  /// collect data for a single walker quantity of complex array type into the current buffer row
  template<unsigned D>
  inline void collect(const std::string& name, Array<std::complex<T>, D> arr)
  {
    ScopedTimer timer(walker_log_buffer_timers_[Timer::COLLECT]);
    size_t n1 = arr.size(0);
    size_t n2, n3, n4;
    n2 = n3 = n4 = 0;
    if (D > 4)
      throw std::runtime_error("WalkerLogBuffer::collect  Only arrays up to dimension 4 are currently supported.");
    if (D > 1)
      n2 = arr.size(1);
    if (D > 2)
      n3 = arr.size(2);
    if (D > 3)
      n4 = arr.size(3);
    size_t irow = 0;
    if (first_collect)
    { // cache walker quantity info on first collect
      WalkerQuantityInfo wqi_(name, 2, walker_data_size, n1, n2, n3, n4);
      quantity_info.push_back(wqi_);
      walker_data_size = wqi_.buffer_end;
      resetRowSize(walker_data_size);
    }
    else
    { // make a new buffer row if needed
      if (quantity_index == 0)
        makeNewRow();
      irow = buffer.size(0) - 1;
    }
    // place the complex array walker quantity into the current buffer row
    auto& wqi   = quantity_info[quantity_index];
    auto& arr1d = arr.storage();
    size_t n    = 0;
    for (size_t i = 0; i < arr1d.size(); ++i)
    {
      buffer(irow, wqi.buffer_start + n) = std::real(arr1d[i]);
      ++n;
      buffer(irow, wqi.buffer_start + n) = std::imag(arr1d[i]);
      ++n;
    }
    quantity_index++;
  }


  /// add a data row from another buffer to this one
  inline void addRow(WalkerLogBuffer<T> other, size_t i)
  {
    ScopedTimer timer(walker_log_buffer_timers_[Timer::ADD_ROW]);
    auto& other_buffer = other.buffer;
    if (first_collect)
    {
      resetRowSize(other_buffer.size(1));
      quantity_info = other.quantity_info;
      first_collect = false;
    }
    else
    {
      if (buffer.size(1) != other_buffer.size(1))
        throw std::runtime_error("WalkerLogBuffer::add_row  Row sizes must match.");
      makeNewRow();
    }
    size_t ib = buffer.size(0) - 1;
    for (size_t j = 0; j < buffer.size(1); ++j)
      buffer(ib, j) = other_buffer(i, j);
  }


  /// write a summary of quantities in the buffer
  inline void writeSummary(std::string pad = "  ")
  {
    ScopedTimer timer(walker_log_buffer_timers_[Timer::WRITE]);

    std::string pad2 = pad + "  ";
    std::string pad3 = pad2 + "  ";
    app_log() << std::endl;
    app_log() << pad << "WalkerLogBuffer(" << label << ")" << std::endl;
    app_log() << pad2 << "nrows       = " << buffer.size(0) << std::endl;
    app_log() << pad2 << "row_size    = " << buffer.size(1) << std::endl;
    for (size_t n = 0; n < quantity_info.size(); ++n)
    {
      auto& wqi = quantity_info[n];
      app_log() << pad2 << "quantity " << n << ":  " << wqi.dimension << "  " << wqi.size << "  " << wqi.unit_size
                << "  " << wqi.buffer_start << "  " << wqi.buffer_end << " (" << wqi.name << ")" << std::endl;
    }
    app_log() << pad << "end WalkerLogBuffer(" << label << ")" << std::endl;
  }

  /// write the data_layout for all walker quantities in the HDF file
  inline void registerHDFData(hdf_archive& f)
  {
    auto& top = label;
    f.push(top);
    f.push("data_layout");
    for (auto& wqi : quantity_info)
    {
      f.push(wqi.name);
      f.write(wqi.dimension, "dimension");
      f.write(wqi.shape, "shape");
      f.write(wqi.size, "size");
      f.write(wqi.unit_size, "unit_size");
      f.write(wqi.buffer_start, "index_start");
      f.write(wqi.buffer_end, "index_end");
      f.pop();
    }
    f.pop();
    f.pop();
    if (!f.open_groups())
      throw std::runtime_error("WalkerLogBuffer(" + label +
                               ")::register_hdf_data() some hdf groups are still open at the end of registration");
    hdf_file_pointer = 0;
  }


  /// write the buffer data into the HDF file
  inline void writeHDF(hdf_archive& f) { writeHDF(f, hdf_file_pointer); }


  /// write the buffer data into the HDF file
  inline void writeHDF(hdf_archive& f, hsize_t& file_pointer)
  {
    ScopedTimer timer(walker_log_buffer_timers_[Timer::WRITE]);
    auto& top = label;
    hsize_t dims[2];
    dims[0] = buffer.size(0);
    dims[1] = buffer.size(1);
    if (dims[0] > 0)
    {
      f.push(top);
      h5d_append(f.top(), "data", file_pointer, buffer.dim(), dims, buffer.data());
      f.pop();
    }
    f.flush();
  }

private:
  /// make space as quantities are added to the buffer for the first time
  inline void resetRowSize(size_t row_size)
  {
    ScopedTimer timer(walker_log_buffer_timers_[Timer::RESET]);
    auto nrows = buffer.size(0);
    if (nrows == 0)
      nrows++;
    if (nrows != 1)
      throw std::runtime_error("WalkerLogBuffer::reset_rowsize  row_size (number of columns) should only be changed "
                               "during growth of the first row.");
    auto buffer_old(buffer);
    buffer.resize(nrows, row_size);
    std::copy_n(buffer_old.data(), buffer_old.size(), buffer.data());
    if (buffer.size(0) != 1)
      throw std::runtime_error(
          "WalkerLogBuffer::reset_rowsize  buffer should contain only a single row upon completion.");
    if (buffer.size(1) != row_size)
      throw std::runtime_error("WalkerLogBuffer::reset_rowsize  length of buffer row should match the requested "
                               "row_size following the reset/udpate.");
  }

  /// allocate a full new row at the end of the buffer
  inline void makeNewRow()
  {
    ScopedTimer timer(walker_log_buffer_timers_[Timer::MAKE_NEW_ROW]);

    size_t nrows    = buffer.size(0);
    size_t row_size = buffer.size(1);
    if (row_size == 0)
      throw std::runtime_error("WalkerLogBuffer::makeNewRow  Cannot make a new row of size zero.");
    nrows++;
    // resizing buffer(type Array) doesn't preserve data. Thus keep old data and copy over
    auto buffer_old(buffer);
    buffer.resize(nrows, row_size);
    std::copy_n(buffer_old.data(), buffer_old.size(), buffer.data());
  }
};

extern template class WalkerLogBuffer<WLog::Int>;
extern template class WalkerLogBuffer<WLog::Real>;

} // namespace qmcplusplus

#endif
