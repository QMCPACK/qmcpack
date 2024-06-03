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

#include "WalkerTraceManager.h"
#include "WalkerTraceInput.h"
#include "Concurrency/OpenMP.h"


namespace qmcplusplus
{


WalkerTraceManager::WalkerTraceManager(WalkerTraceInput& inp, bool allow_traces, std::string series_root, Communicate* comm)
{
  state.reset();
  communicator              = comm;
  file_root                 = series_root;
  bool driver_allows_traces = allow_traces; // driver allows traces or not

  bool traces_requested     = inp.present;  // xml input present or not
  // determine whether walker traces will be active
  state.traces_active       = traces_requested && driver_allows_traces;

  if (state.traces_active)
  {
    if (omp_get_thread_num() == 0)
    {
      app_log() << "\n  WalkerTraceManager::put() " << std::endl;
      app_log() << "    traces requested      : " << traces_requested << std::endl;
      app_log() << "    driver allows traces  : " << driver_allows_traces << std::endl;
      app_log() << "    traces active         : " << state.traces_active << std::endl;
      app_log() << std::endl;
    }
    // retrieve input data
    state.step_period   = inp.get<int>("step_period");
    state.verbose       = inp.get<bool>("verbose");
    bool qtiles         = inp.get<bool>("qtiles");
    write_particle_data = inp.get<bool>("particle");
    write_min_data      = inp.get<bool>("min")    && qtiles;
    write_max_data      = inp.get<bool>("max")    && qtiles;
    write_med_data      = inp.get<bool>("median") && qtiles;
  }
  
  // label min energy walker buffers for HDF file write
  wmin_property_int_buffer.label  = "wmin_property_int";
  wmin_property_real_buffer.label = "wmin_property_real";
  wmin_particle_real_buffer.label = "wmin_particle_real";
  
  // label max energy walker buffers for HDF file write
  wmax_property_int_buffer.label  = "wmax_property_int";
  wmax_property_real_buffer.label = "wmax_property_real";
  wmax_particle_real_buffer.label = "wmax_particle_real";
  
  // label median energy walker buffers for HDF file write
  wmed_property_int_buffer.label  = "wmed_property_int";
  wmed_property_real_buffer.label = "wmed_property_real";
  wmed_particle_real_buffer.label = "wmed_particle_real";

  // walker quantity information ("data_layout") has not be put in the HDF file yet
  registered_hdf = false;
}


WalkerTraceCollector* WalkerTraceManager::makeCollector()
{
  if (state.verbose) app_log() << "WalkerTraceManager::makeCollector " << std::endl;
  WalkerTraceCollector* tc = new WalkerTraceCollector();
  tc->state = state;
  return tc;
}


void WalkerTraceManager::startRun(std::vector<WalkerTraceCollector*>& collectors)
{
  if (!state.traces_active) return; // no-op for driver if traces are inactive
  if (state.verbose) app_log() << "WalkerTraceManager::startRun " << std::endl;
  // transfer step_period, verbosity, etc settings to trace collectors
  for (auto& tc: collectors)
    tc->state = state;
  // check data size consistency among the trace collector buffers
  checkCollectors(collectors);
  // open the traces file
  openFile(collectors);
}


void WalkerTraceManager::stopRun()
{
  if (!state.traces_active) return; // no-op for driver if traces are inactive
  if (state.verbose) app_log() << "WalkerTraceManager::stopRun " << std::endl;
  // close the traces file
  closeFile();
}


void WalkerTraceManager::writeBuffers(std::vector<WalkerTraceCollector*>& collectors)
{
  if (!state.traces_active) return; // no-op for driver if traces are inactive
  if (state.verbose) app_log() << "WalkerTraceManager::writeBuffers "<<std::endl;

  if(write_min_data)
  {// resize min energy walker buffers to zero rows
    wmin_property_int_buffer.resetBuffer();
    wmin_property_real_buffer.resetBuffer();
    wmin_particle_real_buffer.resetBuffer();
  }
  if(write_max_data)
  {// resize max energy walker buffers to zero rows
    wmax_property_int_buffer.resetBuffer();
    wmax_property_real_buffer.resetBuffer();
    wmax_particle_real_buffer.resetBuffer();
  }
  if(write_med_data)
  {// resize median energy walker buffers to zero rows
    wmed_property_int_buffer.resetBuffer();
    wmed_property_real_buffer.resetBuffer();
    wmed_particle_real_buffer.resetBuffer();
  }

  // collect energy information and extract info from min/max/median energy walkers
  if(write_min_data || write_max_data || write_med_data)
  {
    // gather per energy and step data for all walker throughout the MC block
    for (size_t c=0; c<collectors.size(); ++c)
    {
      auto& tc = *collectors[c];
      tc.checkBuffers();
      for (size_t r=0; r<tc.energies.size(); ++r)
        energy_order.push_back(std::make_tuple(tc.steps[r],tc.energies[r],c,r));
    }
    // sort the data by step and energy to enable selection of min/max/median energy walker data
    std::sort(energy_order.begin(), energy_order.end());
    // select out the min/max/median energy walker data and store in rank-level buffers
    size_t n=0;
    size_t n1,n2;
    size_t prev_step;
    for (auto& v: energy_order)
    {
      auto step = std::get<0>(v);
      if (n==0)
      {
        n1 = n;
        prev_step = step;
      }
      if (step!=prev_step || n==energy_order.size()-1)
      {
        // for a given step, find data for min/max/median energy walkers
        //   n1/n2 are indices of the first/last data in energy_order for this step
        if (step!=prev_step) n2 = n-1;
        if (n==energy_order.size()-1) n2 = n;
        auto nmin = n1;        // index of minimum energy walker for this step
        auto nmax = n2;        // index of maximum energy walker for this step
        auto nmed = (n1+n2)/2; // index of median  energy walker for this step
        size_t c,r;
        if(write_min_data)
        {//  cache data for minimum energy walker
          c = std::get<2>(energy_order[nmin]);
          r = std::get<3>(energy_order[nmin]);
          wmin_property_int_buffer.addRow(collectors[c]->walker_property_int_buffer,r);
          wmin_property_real_buffer.addRow(collectors[c]->walker_property_real_buffer,r);
          wmin_particle_real_buffer.addRow(collectors[c]->walker_particle_real_buffer,r);
        }
        if(write_max_data)
        {//  cache data for maximum energy walker
          c = std::get<2>(energy_order[nmax]);
          r = std::get<3>(energy_order[nmax]);
          wmax_property_int_buffer.addRow(collectors[c]->walker_property_int_buffer,r);
          wmax_property_real_buffer.addRow(collectors[c]->walker_property_real_buffer,r);
          wmax_particle_real_buffer.addRow(collectors[c]->walker_particle_real_buffer,r);
        }
        if(write_med_data)
        {//  cache data for median energy walker
          c = std::get<2>(energy_order[nmed]);
          r = std::get<3>(energy_order[nmed]);
          wmed_property_int_buffer.addRow(collectors[c]->walker_property_int_buffer,r);
          wmed_property_real_buffer.addRow(collectors[c]->walker_property_real_buffer,r);
          wmed_particle_real_buffer.addRow(collectors[c]->walker_particle_real_buffer,r);
        }
        // reset pointers
        n1 = n;
        prev_step = step;
      }
      n++;
    }
    energy_order.resize(0);
  }

  // write buffer data to file
  writeBuffersHDF(collectors);
}


void WalkerTraceManager::checkCollectors(std::vector<WalkerTraceCollector*>& collectors)
{
  if (state.verbose) app_log() << "WalkerTraceManager::checkCollectors" << std::endl;
  if (collectors.size() > 0)
  {
    bool all_same = true;
    WalkerTraceCollector& ref = *collectors[0];
    for (int i = 0; i < collectors.size(); ++i)
    {
      WalkerTraceCollector& tc = *collectors[i];
      all_same &= tc.walker_property_int_buffer.sameAs(ref.walker_property_int_buffer);
      all_same &= tc.walker_property_real_buffer.sameAs(ref.walker_property_real_buffer);
      all_same &= tc.walker_particle_real_buffer.sameAs(ref.walker_particle_real_buffer);
    }
    if (!all_same)
    {
      throw std::runtime_error("WalkerTraceManager::checkCollectors  trace buffer widths of collectors do not match\n  contiguous write is "
                               "impossible\n  this was first caused by collectors contributing array traces from identical, but "
                               "differently named, particlesets such as e, e2, e3 ... (fixed)\n  please check the WalkerTraceManager "
                               "summaries printed above");
    }
  }
}


void WalkerTraceManager::openFile(std::vector<WalkerTraceCollector*>& collectors)
{
  if (state.verbose) app_log() << "WalkerTraceManager::openFile "<<std::endl;
  openHDFFile(collectors);
}


void WalkerTraceManager::closeFile()
{
  if (state.verbose) app_log() << "WalkerTraceManager::closeFile " << std::endl;
  closeHDFFile();
}


void WalkerTraceManager::openHDFFile(std::vector<WalkerTraceCollector*>& collectors)
{
  if (state.verbose) app_log() << "WalkerTraceManager::openHDFFile " << std::endl;
  if (collectors.size() == 0) 
    throw std::runtime_error("WalkerTraceManager::openHDFFile  no trace collectors exist, cannot open file");
  // each rank opens a wtraces.h5 file
  int nprocs = communicator->size();
  int rank   = communicator->rank();
  std::array<char, 32> ptoken;
  std::string file_name = file_root;
  if (nprocs > 1)
  {// extend the file name to include processor/rank information
    int length{0};
    if (nprocs > 10000)
      length = std::snprintf(ptoken.data(), ptoken.size(), ".p%05d", rank);
    else if (nprocs > 1000)
      length = std::snprintf(ptoken.data(), ptoken.size(), ".p%04d", rank);
    else
      length = std::snprintf(ptoken.data(), ptoken.size(), ".p%03d", rank);
    if (length < 0)
      throw std::runtime_error("Error generating filename");
    file_name.append(ptoken.data(), length);
  }
  file_name += ".wtraces.h5";
  if (state.verbose) app_log() << "WalkerTraceManager::openHDFFile  opening traces hdf file " << file_name << std::endl;
  // create the hdf archive
  hdf_file        = std::make_unique<hdf_archive>();
  // open the file
  bool successful = hdf_file->create(file_name);
  if (!successful)
    throw std::runtime_error("WalkerTraceManager::openHDFFile  failed to open hdf file " + file_name);
}


void WalkerTraceManager::writeBuffersHDF(std::vector<WalkerTraceCollector*>& collectors)
{
  if (state.verbose) app_log() << "WalkerTraceManager::writeBuffersHDF " << std::endl;
  WalkerTraceCollector& tc_lead = *collectors[0];
  if(!registered_hdf)
  {// write walker quantity information ("data_layout") for each buffer in the HDF file
    //  create data_layout for all-walker buffers
    tc_lead.walker_property_int_buffer.registerHDFData(*hdf_file);
    tc_lead.walker_property_real_buffer.registerHDFData(*hdf_file);
    if(write_particle_data)
      tc_lead.walker_particle_real_buffer.registerHDFData(*hdf_file);
    if(write_min_data)
    {//  create data_layout for min energy walker buffers
      wmin_property_int_buffer.registerHDFData(*hdf_file);
      wmin_property_real_buffer.registerHDFData(*hdf_file);
      wmin_particle_real_buffer.registerHDFData(*hdf_file);
    }
    if(write_max_data)
    {//  create data_layout for max energy walker buffers
      wmax_property_int_buffer.registerHDFData(*hdf_file);
      wmax_property_real_buffer.registerHDFData(*hdf_file);
      wmax_particle_real_buffer.registerHDFData(*hdf_file);
    }
    if(write_med_data)
    {//  create data_layout for median energy walker buffers
      wmed_property_int_buffer.registerHDFData(*hdf_file);
      wmed_property_real_buffer.registerHDFData(*hdf_file);
      wmed_particle_real_buffer.registerHDFData(*hdf_file);
    }
    // walker quantity information ("data_layout") has now been added to HDF, do not repeat
    registered_hdf = true;
  }
  // write data for all-walker buffers to HDF
  for (int ip = 0; ip < collectors.size(); ++ip)
  {
    WalkerTraceCollector& tc = *collectors[ip];
    tc.walker_property_int_buffer.writeHDF(*hdf_file, tc_lead.walker_property_int_buffer.hdf_file_pointer);
    tc.walker_property_real_buffer.writeHDF(*hdf_file, tc_lead.walker_property_real_buffer.hdf_file_pointer);
    if(write_particle_data)
      tc.walker_particle_real_buffer.writeHDF(*hdf_file, tc_lead.walker_particle_real_buffer.hdf_file_pointer);
  }
  if(write_min_data)
  {// write data for min energy walker buffers to HDF
    wmin_property_int_buffer.writeHDF(*hdf_file);
    wmin_property_real_buffer.writeHDF(*hdf_file);
    wmin_particle_real_buffer.writeHDF(*hdf_file);
  }
  if(write_max_data)
  {// write data for max energy walker buffers to HDF
    wmax_property_int_buffer.writeHDF(*hdf_file);
    wmax_property_real_buffer.writeHDF(*hdf_file);
    wmax_particle_real_buffer.writeHDF(*hdf_file);
  }
  if(write_med_data)
  {// write data for median energy walker buffers to HDF
    wmed_property_int_buffer.writeHDF(*hdf_file);
    wmed_property_real_buffer.writeHDF(*hdf_file);
    wmed_particle_real_buffer.writeHDF(*hdf_file);
  }
}


void WalkerTraceManager::closeHDFFile()
{ 
  if (state.verbose) app_log() << "WalkerTraceManager::closeHDFFile " << std::endl;
  hdf_file.reset(); 
}



} // namespace qmcplusplus
