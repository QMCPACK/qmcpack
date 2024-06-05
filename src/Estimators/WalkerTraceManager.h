//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_WALKERTRACEMANAGER_H
#define QMCPLUSPLUS_WALKERTRACEMANAGER_H

#include "WalkerTraceCollector.h"


namespace qmcplusplus
{


struct WalkerTraceInput;


/** Driver-level resource for walker trace collection.
 *
 *    This class manages the HDF file (open/close/write).
 *    Walker buffer data from all crowd-level WalkerTraceCollectors are written
 *    to the HDF file at the end of each MC block.
 *
 *    Just prior to the write, this class examines the distribution of walker
 *    energies on its rank for each MC step in the MC block, identifies the 
 *    minimum/maximum/median energy walkers and buffers their full data 
 *    (including per-particle information) for the write.
 *    This data corresponds to walker information at specific quantiles of the 
 *    energy distribution.
 *    See WalkerTraceManager::writeBuffers()
 *
 *    Writing per-particle information for all walkers is optional, as is 
 *    writing data for the minimum/maximum/median energy walkers.
 *    Scalar "property" data for each walker is always written.
 */
class WalkerTraceManager
{
public:
  WalkerTraceState state;

private:
  /// file prefix for the current driver
  std::string file_root;
  Communicate* communicator;
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
  std::vector<std::tuple<size_t,WTrace::Real,size_t,size_t>> energy_order;

  /// buffer containing integer properties for the minimum energy walkers
  WalkerTraceBuffer<WTrace::Int>  wmin_property_int_buffer;
  /// buffer containing real-valued properties for the minimum energy walkers
  WalkerTraceBuffer<WTrace::Real> wmin_property_real_buffer;
  /// buffer containing per-particle properties for the minimum energy walkers
  WalkerTraceBuffer<WTrace::Real> wmin_particle_real_buffer;

  /// buffer containing integer properties for the maximum energy walkers
  WalkerTraceBuffer<WTrace::Int>  wmax_property_int_buffer;
  /// buffer containing real-valued properties for the maximum energy walkers
  WalkerTraceBuffer<WTrace::Real> wmax_property_real_buffer;
  /// buffer containing per-particle properties for the maximum energy walkers
  WalkerTraceBuffer<WTrace::Real> wmax_particle_real_buffer;

  /// buffer containing integer properties for the median energy walkers
  WalkerTraceBuffer<WTrace::Int>  wmed_property_int_buffer;
  /// buffer containing real-valued properties for the median energy walkers
  WalkerTraceBuffer<WTrace::Real> wmed_property_real_buffer;
  /// buffer containing per-particle properties for the median energy walkers
  WalkerTraceBuffer<WTrace::Real> wmed_particle_real_buffer;

public:
  WalkerTraceManager(WalkerTraceInput& inp, bool allow_traces, std::string series_root, Communicate* comm = 0);

  /// create a WalkerTraceCollector (legacy drivers only, "cloning" style)
  WalkerTraceCollector* makeCollector();

  /// open the traces file and check consistency of the collectors at the start of a run
  void startRun(std::vector<WalkerTraceCollector*>& collectors);

  /// close the traces file at the end of a run
  void stopRun();

  /// collect min/max/median walker data and write buffered walker trace data to file
  void writeBuffers(std::vector<WalkerTraceCollector*>& collectors);

private:
  /// check consistency of walker buffer row sizes
  void checkCollectors(std::vector<WalkerTraceCollector*>& collectors);

  /// open the traces file
  void openFile(std::vector<WalkerTraceCollector*>& collectors);

  /// close the traces file
  void closeFile();

  /// open the traces file (HDF format)
  void openHDFFile(std::vector<WalkerTraceCollector*>& collectors);

  /// write data buffers to the traces file (HDF format)
  void writeBuffersHDF(std::vector<WalkerTraceCollector*>& collectors);

  /// close the traces file (HDF format)
  void closeHDFFile();
};



} // namespace qmcplusplus

#endif
