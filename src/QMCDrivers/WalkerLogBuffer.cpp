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

template class WalkerLogBuffer<WLog::Int>;
template class WalkerLogBuffer<WLog::Real>;

} // namespace qmcplusplus
