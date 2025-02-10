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

template class WalkerLogBuffer<WLog::Int>;
template class WalkerLogBuffer<WLog::Real>;

template void WalkerLogBuffer<WLog::Int>::collect(const std::string& name, Array<WLog::Int, 2>);
template void WalkerLogBuffer<WLog::Real>::collect<2>(const std::string& name, Array<WLog::Real, 2>);
template void WalkerLogBuffer<WLog::Real>::collect<1>(const std::string& name, Array<WLog::Real, 1>);
template void WalkerLogBuffer<WLog::Real>::collect<2>(const std::string& name,
                                                             Array<std::complex<WLog::Real>, 2>);
template void WalkerLogBuffer<WLog::Real>::collect<1>(const std::string& name,
                                                             Array<std::complex<WLog::Real>, 1>);

} // namespace qmcplusplus
