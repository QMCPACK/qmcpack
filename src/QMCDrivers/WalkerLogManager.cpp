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

#include "WalkerLogManager.h"
#include "WalkerLogInput.h"
#include "Concurrency/Info.hpp"


namespace qmcplusplus
{

TimerNameList_t<WalkerLogManager::Timer> WalkerLogManager::create_names(const std::string_view& my_name)
{
  TimerNameList_t<Timer> timer_names;
  using namespace std::string_literals;
  std::string prefix{"WalkerLog:"s + std::string{my_name} + "::"s};
  for (std::size_t i = 0; i < suffixes_.size(); ++i)
    timer_names.push_back({static_cast<Timer>(i), prefix + std::string{suffixes_[i]}});
  return timer_names;
}

WalkerLogManager::WalkerLogManager(WalkerLogInput& inp, bool allow_logs, std::string series_root, Communicate* comm)
    : walker_log_timers_(getGlobalTimerManager(), create_names(my_name_), timer_level_medium)
{
  communicator            = comm;
  file_root               = series_root;
  bool driver_allows_logs = allow_logs; // driver allows logs or not

  bool logs_requested = inp.present; // xml input present or not
  // determine whether walker logs will be active
  state.logs_active = logs_requested && driver_allows_logs;

  if (state.logs_active)
  {
    if (Concurrency::getWorkerId() == 0)
    {
      app_log() << "\n  WalkerLogManager::put() " << std::endl;
      app_log() << "    logs requested      : " << logs_requested << std::endl;
      app_log() << "    driver allows logs  : " << driver_allows_logs << std::endl;
      app_log() << "    logs active         : " << state.logs_active << std::endl;
      app_log() << std::endl;
    }
    // retrieve input data
    state.step_period   = inp.get<int>("step_period");
    state.verbose       = inp.get<bool>("verbose");
    bool quantiles      = inp.get<bool>("quantiles");
    write_particle_data = inp.get<bool>("particle");
    write_min_data      = inp.get<bool>("min") && quantiles;
    write_max_data      = inp.get<bool>("max") && quantiles;
    write_med_data      = inp.get<bool>("median") && quantiles;
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


std::unique_ptr<WalkerLogCollector> WalkerLogManager::makeCollector() const
{
  if (state.verbose)
    app_log() << "WalkerLogManager::makeCollector " << std::endl;
  return std::make_unique<WalkerLogCollector>(state);
}


void WalkerLogManager::startRun(RefVector<WalkerLogCollector>&& collectors)
{
  if (!state.logs_active)
    return; // no-op for driver if logs are inactive

  ScopedTimer timer(walker_log_timers_[Timer::START]);

  collectors_in_run_ = std::move(collectors);
  if (collectors_in_run_.empty())
    throw std::runtime_error("BUG collectors are empty but walker logs are active");
  if (state.verbose)
    app_log() << "WalkerLogManager::startRun " << std::endl;
  // check data size consistency among the log collector buffers
  checkCollectors(collectors_in_run_);
  // open the logs file
  openFile(collectors_in_run_);
}


void WalkerLogManager::stopRun()
{
  if (!state.logs_active)
    return; // no-op for driver if logs are inactive

  ScopedTimer timer(walker_log_timers_[Timer::STOP]);

  if (state.verbose)
    app_log() << "WalkerLogManager::stopRun " << std::endl;
  collectors_in_run_.clear();
  // close the logs file
  closeFile();
}


void WalkerLogManager::writeBuffers()
{
  if (!state.logs_active)
    return; // no-op for driver if logs are inactive

  ScopedTimer timer(walker_log_timers_[Timer::WRITE]);

  if (state.verbose)
    app_log() << "WalkerLogManager::writeBuffers " << std::endl;

  if (write_min_data)
  { // resize min energy walker buffers to zero rows
    wmin_property_int_buffer.resetBuffer();
    wmin_property_real_buffer.resetBuffer();
    wmin_particle_real_buffer.resetBuffer();
  }
  if (write_max_data)
  { // resize max energy walker buffers to zero rows
    wmax_property_int_buffer.resetBuffer();
    wmax_property_real_buffer.resetBuffer();
    wmax_particle_real_buffer.resetBuffer();
  }
  if (write_med_data)
  { // resize median energy walker buffers to zero rows
    wmed_property_int_buffer.resetBuffer();
    wmed_property_real_buffer.resetBuffer();
    wmed_particle_real_buffer.resetBuffer();
  }

  const RefVector<WalkerLogCollector>& collectors = collectors_in_run_;
  // collect energy information and extract info from min/max/median energy walkers
  if (write_min_data || write_max_data || write_med_data)
  {
    // gather per energy and step data for all walker throughout the MC block
    for (size_t c = 0; c < collectors.size(); ++c)
    {
      WalkerLogCollector& tc = collectors[c];
      tc.checkBuffers();
      for (size_t r = 0; r < tc.energies.size(); ++r)
        energy_order.push_back(std::make_tuple(tc.steps[r], tc.energies[r], c, r));
    }
    // sort the data by step and energy to enable selection of min/max/median energy walker data
    std::sort(energy_order.begin(), energy_order.end());
    // select out the min/max/median energy walker data and store in rank-level buffers
    size_t n = 0;
    size_t n1, n2;
    size_t prev_step;
    for (auto& v : energy_order)
    {
      auto step = std::get<0>(v);
      if (n == 0)
      {
        n1        = n;
        prev_step = step;
      }
      if (step != prev_step || n == energy_order.size() - 1)
      {
        // for a given step, find data for min/max/median energy walkers
        //   n1/n2 are indices of the first/last data in energy_order for this step
        if (step != prev_step)
          n2 = n - 1;
        if (n == energy_order.size() - 1)
          n2 = n;
        auto nmin = n1;            // index of minimum energy walker for this step
        auto nmax = n2;            // index of maximum energy walker for this step
        auto nmed = (n1 + n2) / 2; // index of median  energy walker for this step
        size_t c, r;
        if (write_min_data)
        { //  cache data for minimum energy walker
          c                      = std::get<2>(energy_order[nmin]);
          r                      = std::get<3>(energy_order[nmin]);
          WalkerLogCollector& tc = collectors[c];
          wmin_property_int_buffer.addRow(tc.walker_property_int_buffer, r);
          wmin_property_real_buffer.addRow(tc.walker_property_real_buffer, r);
          wmin_particle_real_buffer.addRow(tc.walker_particle_real_buffer, r);
        }
        if (write_max_data)
        { //  cache data for maximum energy walker
          c                      = std::get<2>(energy_order[nmax]);
          r                      = std::get<3>(energy_order[nmax]);
          WalkerLogCollector& tc = collectors[c];
          wmax_property_int_buffer.addRow(tc.walker_property_int_buffer, r);
          wmax_property_real_buffer.addRow(tc.walker_property_real_buffer, r);
          wmax_particle_real_buffer.addRow(tc.walker_particle_real_buffer, r);
        }
        if (write_med_data)
        { //  cache data for median energy walker
          c                      = std::get<2>(energy_order[nmed]);
          r                      = std::get<3>(energy_order[nmed]);
          WalkerLogCollector& tc = collectors[c];
          wmed_property_int_buffer.addRow(tc.walker_property_int_buffer, r);
          wmed_property_real_buffer.addRow(tc.walker_property_real_buffer, r);
          wmed_particle_real_buffer.addRow(tc.walker_particle_real_buffer, r);
        }
        // reset pointers
        n1        = n;
        prev_step = step;
      }
      n++;
    }
    energy_order.resize(0);
  }

  // write buffer data to file
  writeBuffersHDF();
}


void WalkerLogManager::checkCollectors(const RefVector<WalkerLogCollector>& collectors) const
{
  if (state.verbose)
    app_log() << "WalkerLogManager::checkCollectors" << std::endl;
  if (collectors.size() > 0)
  {
    bool all_same           = true;
    WalkerLogCollector& ref = collectors[0];
    for (int i = 0; i < collectors.size(); ++i)
    {
      WalkerLogCollector& tc = collectors[i];
      all_same &= tc.walker_property_int_buffer.sameAs(ref.walker_property_int_buffer);
      all_same &= tc.walker_property_real_buffer.sameAs(ref.walker_property_real_buffer);
      all_same &= tc.walker_particle_real_buffer.sameAs(ref.walker_particle_real_buffer);
    }
    if (!all_same)
    {
      throw std::runtime_error(
          "WalkerLogManager::checkCollectors  log buffer widths of collectors do not match\n  contiguous write is "
          "impossible\n  this was first caused by collectors contributing array logs from identical, but "
          "differently named, particlesets such as e, e2, e3 ... (fixed)\n  please check the WalkerLogManager "
          "summaries printed above");
    }
  }
}


void WalkerLogManager::openFile(const RefVector<WalkerLogCollector>& collectors)
{
  if (state.verbose)
    app_log() << "WalkerLogManager::openFile " << std::endl;
  openHDFFile(collectors);
}


void WalkerLogManager::closeFile()
{
  if (state.verbose)
    app_log() << "WalkerLogManager::closeFile " << std::endl;
  closeHDFFile();
}


void WalkerLogManager::openHDFFile(const RefVector<WalkerLogCollector>& collectors)
{
  if (state.verbose)
    app_log() << "WalkerLogManager::openHDFFile " << std::endl;
  if (collectors.size() == 0)
    throw std::runtime_error("WalkerLogManager::openHDFFile  no log collectors exist, cannot open file");
  // each rank opens a wlogs.h5 file
  int nprocs = communicator->size();
  int rank   = communicator->rank();
  std::array<char, 32> ptoken;
  std::string file_name = file_root;
  if (nprocs > 1)
  { // extend the file name to include processor/rank information
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
  file_name += ".wlogs.h5";
  if (state.verbose)
    app_log() << "WalkerLogManager::openHDFFile  opening logs hdf file " << file_name << std::endl;
  // create the hdf archive
  hdf_file = std::make_unique<hdf_archive>();
  // open the file
  bool successful = hdf_file->create(file_name);
  if (!successful)
    throw std::runtime_error("WalkerLogManager::openHDFFile  failed to open hdf file " + file_name);
}


void WalkerLogManager::writeBuffersHDF()
{
  const RefVector<WalkerLogCollector>& collectors = collectors_in_run_;
  if (state.verbose)
    app_log() << "WalkerLogManager::writeBuffersHDF " << std::endl;
  WalkerLogCollector& tc_lead = collectors[0];
  if (!registered_hdf)
  { // write walker quantity information ("data_layout") for each buffer in the HDF file
    //  create data_layout for all-walker buffers
    tc_lead.walker_property_int_buffer.registerHDFData(*hdf_file);
    tc_lead.walker_property_real_buffer.registerHDFData(*hdf_file);
    if (write_particle_data)
      tc_lead.walker_particle_real_buffer.registerHDFData(*hdf_file);
    if (write_min_data)
    { //  create data_layout for min energy walker buffers
      wmin_property_int_buffer.registerHDFData(*hdf_file);
      wmin_property_real_buffer.registerHDFData(*hdf_file);
      wmin_particle_real_buffer.registerHDFData(*hdf_file);
    }
    if (write_max_data)
    { //  create data_layout for max energy walker buffers
      wmax_property_int_buffer.registerHDFData(*hdf_file);
      wmax_property_real_buffer.registerHDFData(*hdf_file);
      wmax_particle_real_buffer.registerHDFData(*hdf_file);
    }
    if (write_med_data)
    { //  create data_layout for median energy walker buffers
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
    WalkerLogCollector& tc = collectors[ip];
    tc.walker_property_int_buffer.writeHDF(*hdf_file, tc_lead.walker_property_int_buffer.hdf_file_pointer);
    tc.walker_property_real_buffer.writeHDF(*hdf_file, tc_lead.walker_property_real_buffer.hdf_file_pointer);
    if (write_particle_data)
      tc.walker_particle_real_buffer.writeHDF(*hdf_file, tc_lead.walker_particle_real_buffer.hdf_file_pointer);
  }
  if (write_min_data)
  { // write data for min energy walker buffers to HDF
    wmin_property_int_buffer.writeHDF(*hdf_file);
    wmin_property_real_buffer.writeHDF(*hdf_file);
    wmin_particle_real_buffer.writeHDF(*hdf_file);
  }
  if (write_max_data)
  { // write data for max energy walker buffers to HDF
    wmax_property_int_buffer.writeHDF(*hdf_file);
    wmax_property_real_buffer.writeHDF(*hdf_file);
    wmax_particle_real_buffer.writeHDF(*hdf_file);
  }
  if (write_med_data)
  { // write data for median energy walker buffers to HDF
    wmed_property_int_buffer.writeHDF(*hdf_file);
    wmed_property_real_buffer.writeHDF(*hdf_file);
    wmed_particle_real_buffer.writeHDF(*hdf_file);
  }
}


void WalkerLogManager::closeHDFFile()
{
  if (state.verbose)
    app_log() << "WalkerLogManager::closeHDFFile " << std::endl;
  hdf_file.reset();
}


} // namespace qmcplusplus
