//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "WalkerLogBuffer.h"

namespace qmcplusplus
{

template<typename T>
WalkerLogBuffer<T>::WalkerLogBuffer()
    : walker_log_buffer_timers_(getGlobalTimerManager(), create_names(my_name_), timer_level_medium)
{
  label            = "?";
  first_collect    = true;
  walker_data_size = 0;
  quantity_index   = 0;
  resetBuffer();
}

template<typename T>
void WalkerLogBuffer<T>::addRow(WalkerLogBuffer<T> other, size_t i)
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

/// collect data for a single walker quantity of scalar type into the current buffer row
template<typename T>
void WalkerLogBuffer<T>::collect(const std::string& name, const T& value)
{
  ScopedTimer timer(walker_log_buffer_timers_[Timer::COLLECT]);

  size_t irow = 0;
  if (first_collect)
  { // cache walker quantity info on first collect
    WalkerQuantityInfo wqi_(name, 1, walker_data_size);
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
  // place the scalar walker quantity into the current buffer row
  auto& wqi                      = quantity_info[quantity_index];
  buffer(irow, wqi.buffer_start) = value;
  quantity_index++;
}

template<typename T>
template<unsigned D>
inline void WalkerLogBuffer<T>::collect(const std::string& name, Array<T, D> arr)
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

template<typename T>
template<unsigned D>
inline void WalkerLogBuffer<T>::collect(const std::string& name, Array<std::complex<T>, D> arr)
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

template<typename T>
void WalkerLogBuffer<T>::makeNewRow()
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

template<typename T>
void WalkerLogBuffer<T>::registerHDFData(hdf_archive& f)
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

template<typename T>
void WalkerLogBuffer<T>::resetRowSize(size_t row_size)
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

template<typename T>
void WalkerLogBuffer<T>::writeHDF(hdf_archive& f)
{
  writeHDF(f, hdf_file_pointer);
}

template<typename T>
void WalkerLogBuffer<T>::writeHDF(hdf_archive& f, hsize_t& file_pointer)
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

template<typename T>
void WalkerLogBuffer<T>::writeSummary(std::string pad)
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
    app_log() << pad2 << "quantity " << n << ":  " << wqi.dimension << "  " << wqi.size << "  " << wqi.unit_size << "  "
              << wqi.buffer_start << "  " << wqi.buffer_end << " (" << wqi.name << ")" << std::endl;
  }
  app_log() << pad << "end WalkerLogBuffer(" << label << ")" << std::endl;
}

template class WalkerLogBuffer<WLog::Int>;
template class WalkerLogBuffer<WLog::Real>;

template void WalkerLogBuffer<WLog::Int>::collect(const std::string& name, Array<WLog::Int, 2>);
template void WalkerLogBuffer<WLog::Real>::collect<2>(const std::string& name, Array<WLog::Real, 2>);
template void WalkerLogBuffer<WLog::Real>::collect<1>(const std::string& name, Array<WLog::Real, 1>);
template void WalkerLogBuffer<WLog::Real>::collect<2>(const std::string& name, Array<std::complex<WLog::Real>, 2>);
template void WalkerLogBuffer<WLog::Real>::collect<1>(const std::string& name, Array<std::complex<WLog::Real>, 1>);

} // namespace qmcplusplus
