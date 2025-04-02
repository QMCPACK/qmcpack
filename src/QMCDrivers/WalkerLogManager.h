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


#ifndef QMCPLUSPLUS_WALKERLOGMANAGER_H
#define QMCPLUSPLUS_WALKERLOGMANAGER_H

#include "WalkerLogCollector.h"
#include "type_traits/template_types.hpp"
#include "Utilities/TimerManager.h"

namespace qmcplusplus
{


struct WalkerLogInput;


/** Driver-level resource for walker log collection.
 *
 *    This class manages the HDF file (open/close/write).
 *    Walker buffer data from all crowd-level WalkerLogCollectors are written
 *    to the HDF file at the end of each MC block.
 *
 *    Just prior to the write, this class examines the distribution of walker
 *    energies on its rank for each MC step in the MC block, identifies the 
 *    minimum/maximum/median energy walkers and buffers their full data 
 *    (including per-particle information) for the write.
 *    This data corresponds to walker information at specific quantiles of the 
 *    energy distribution.
 *    See WalkerLogManager::writeBuffers()
 *
 *    Writing per-particle information for all walkers is optional, as is 
 *    writing data for the minimum/maximum/median energy walkers.
 *    Scalar "property" data for each walker is always written.
 */
class WalkerLogManager
{
private:
  /// file prefix for the current driver
  std::string file_root;
  Communicate* communicator;
  /// output state
  WalkerLogState state;
  /// access to HDF file
  std::unique_ptr<hdf_archive> hdf_file;
  /// whether walker quantity data ("data_layout") has been recorded in HDF
  bool registered_hdf;

  /// whether to write per-particle data for all walkers
  bool write_particle_data;
  /// whether to write full data for the minimum energy walker at each step
  bool write_min_data;
  /// whether to write full data for the maximum energy walker at each step
  bool write_max_data;
  /// whether to write full data for the median energy walker at each step
  bool write_med_data;

  /// used to sort energy information for identifying walkers by energy quantile
  std::vector<std::tuple<size_t, WLog::Real, size_t, size_t>> energy_order;

  /// buffer containing integer properties for the minimum energy walkers
  WalkerLogBuffer<WLog::Int> wmin_property_int_buffer;
  /// buffer containing real-valued properties for the minimum energy walkers
  WalkerLogBuffer<WLog::Real> wmin_property_real_buffer;
  /// buffer containing per-particle properties for the minimum energy walkers
  WalkerLogBuffer<WLog::Real> wmin_particle_real_buffer;

  /// buffer containing integer properties for the maximum energy walkers
  WalkerLogBuffer<WLog::Int> wmax_property_int_buffer;
  /// buffer containing real-valued properties for the maximum energy walkers
  WalkerLogBuffer<WLog::Real> wmax_property_real_buffer;
  /// buffer containing per-particle properties for the maximum energy walkers
  WalkerLogBuffer<WLog::Real> wmax_particle_real_buffer;

  /// buffer containing integer properties for the median energy walkers
  WalkerLogBuffer<WLog::Int> wmed_property_int_buffer;
  /// buffer containing real-valued properties for the median energy walkers
  WalkerLogBuffer<WLog::Real> wmed_property_real_buffer;
  /// buffer containing per-particle properties for the median energy walkers
  WalkerLogBuffer<WLog::Real> wmed_particle_real_buffer;

  RefVector<WalkerLogCollector> collectors_in_run_;

  /// medium granularity timers
  TimerList_t walker_log_timers_;
  static constexpr std::string_view my_name_{"WalkerLogManager"};
  enum Timer
  {
    START = 0,
    STOP,
    WRITE,
  };
  static constexpr std::array<std::string_view, 3> suffixes_{"start", "stop", "write"};
  static TimerNameList_t<Timer> create_names(const std::string_view& my_name);

public:
  WalkerLogManager(WalkerLogInput& inp, bool allow_logs, std::string series_root, Communicate* comm = 0);

  /// create a WalkerLogCollector
  std::unique_ptr<WalkerLogCollector> makeCollector() const;

  /// open the logs file and check consistency of the collectors at the start of a run
  void startRun(RefVector<WalkerLogCollector>&& collectors);

  /// close the logs file at the end of a run
  void stopRun();

  /// collect min/max/median walker data and write buffered walker log data to file
  void writeBuffers();

private:
  /// check consistency of walker buffer row sizes
  void checkCollectors(const RefVector<WalkerLogCollector>& collectors) const;

  /// open the logs file
  void openFile(const RefVector<WalkerLogCollector>& collectors);

  /// close the logs file
  void closeFile();

  /// open the logs file (HDF format)
  void openHDFFile(const RefVector<WalkerLogCollector>& collectors);

  /// write data buffers to the logs file (HDF format)
  void writeBuffersHDF();

  /// close the logs file (HDF format)
  void closeHDFFile();
};


} // namespace qmcplusplus

#endif
