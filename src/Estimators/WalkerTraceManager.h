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

#include <Configuration.h>
#include "WalkerTraceCollector.h"


namespace qmcplusplus
{


class WalkerTraceInput;


class WalkerTraceManager
{
public:
  WalkerTraceState state;

private:
  std::string file_root;
  Communicate* communicator;
  std::unique_ptr<hdf_archive> hdf_file;
  bool registered_hdf;

  // new walker buffers
  bool write_particle_data;
  bool write_min_data;
  bool write_max_data;
  bool write_med_data;

  std::vector<std::tuple<size_t,WTraceReal,size_t,size_t>> energy_order;

  WalkerTraceBuffer<WTraceInt>  wmin_property_int_buffer;
  WalkerTraceBuffer<WTraceReal> wmin_property_real_buffer;
  WalkerTraceBuffer<WTraceReal> wmin_particle_real_buffer;

  WalkerTraceBuffer<WTraceInt>  wmax_property_int_buffer;
  WalkerTraceBuffer<WTraceReal> wmax_property_real_buffer;
  WalkerTraceBuffer<WTraceReal> wmax_particle_real_buffer;

  WalkerTraceBuffer<WTraceInt>  wmed_property_int_buffer;
  WalkerTraceBuffer<WTraceReal> wmed_property_real_buffer;
  WalkerTraceBuffer<WTraceReal> wmed_particle_real_buffer;

public:
  WalkerTraceManager(WalkerTraceInput& inp, bool allow_traces, std::string series_root, Communicate* comm = 0);

  WalkerTraceCollector* makeCollector();

  void startRun(std::vector<WalkerTraceCollector*>& collectors);

  void stopRun();

  //write buffered trace data to file
  void writeBuffers(std::vector<WalkerTraceCollector*>& collectors);

private:
  void checkCollectors(std::vector<WalkerTraceCollector*>& collectors);

  void openFile(std::vector<WalkerTraceCollector*>& collectors);

  void closeFile();

  //hdf file operations
  void openHDFFile(std::vector<WalkerTraceCollector*>& collectors);

  void writeBuffersHDF(std::vector<WalkerTraceCollector*>& collectors);

  void closeHDFFile();
};



} // namespace qmcplusplus

#endif
