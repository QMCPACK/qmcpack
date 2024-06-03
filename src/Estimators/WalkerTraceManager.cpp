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
#include "WalkerTraceManagerInput.h"

#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsData/AttributeSet.h"

#include "Particle/Walker.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"



namespace qmcplusplus
{

using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;


//////////////////////////////////
// WalkerTraceCollector methods //
//////////////////////////////////

WalkerTraceCollector::WalkerTraceCollector()
  : properties_include{"R2Accepted","R2Proposed","LocalEnergy","LocalPotential","Kinetic","ElecElec","ElecIon","LocalECP","NonLocalECP"}
{
  state.reset();
  energy_index = -1;
  steps.resize(0);
  energies.resize(0);
  walker_property_int_buffer.label  = "walker_property_int";
  walker_property_real_buffer.label = "walker_property_real";
  walker_particle_real_buffer.label = "walker_particle_real";
}


void WalkerTraceCollector::startBlock()
{
  if(!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceCollector::startBlock " << std::endl;
  resetBuffers();
}


void WalkerTraceCollector::collect(const MCPWalker& walker, const ParticleSet& pset, const TrialWaveFunction& wfn, const QMCHamiltonian& ham, int step)
{
  if(!state.traces_active) return;

  int current_step = (step==-1) ? pset.current_step : step;
  if(current_step%state.step_period!=0) return;

  auto& bsi = walker_property_int_buffer;
  auto& bsr = walker_property_real_buffer;
  auto& bar = walker_particle_real_buffer;
  
  // collect per-particle walker quantities
  size_t nparticles = walker.R.size();
  size_t ndim       = walker.R[0].size();
  //   per-particle positions (walker.R)
  Rtmp.resize(nparticles,ndim);
  for(size_t p=0;p<nparticles;++p)
    for(size_t d=0;d<ndim;++d)
      Rtmp(p,d) = (WTraceReal)walker.R[p][d];
  bar.collect("R", Rtmp);
  //   per-particle "spin" (walker.spins)
  if (pset.isSpinor())
  {
    Stmp.resize(nparticles);
    for(size_t p=0;p<nparticles;++p)
      Stmp(p) = (WTraceReal)walker.spins[p];
    bar.collect("S", Stmp);
  }
  //   per-particle gradient(log(psi)) (pset.G)
  Gtmp.resize(nparticles,ndim);
  for(size_t p=0;p<nparticles;++p)
    for(size_t d=0;d<ndim;++d)
      Gtmp(p,d) = (WTracePsiVal)pset.G[p][d];
  bar.collect("G", Gtmp);
  //   per-particle laplacian(log(psi)) (pset.L)
  Ltmp.resize(nparticles);
  for(size_t p=0;p<nparticles;++p)
    Ltmp(p) = (WTracePsiVal)pset.L[p];
  bar.collect("L", Ltmp);
  bar.resetCollect();

  // collect integer walker properties
  bsi.collect("step"        , (WTraceInt)current_step         );
  bsi.collect("id"          , (WTraceInt)walker.getWalkerID() );
  bsi.collect("parent_id"   , (WTraceInt)walker.getParentID() );
  bsi.collect("age"         , (WTraceInt)walker.Age           );
  bsi.resetCollect();

  // collect real walker properties
  bsr.collect("weight"      , (WTraceReal)walker.Weight        );
  bsr.collect("multiplicity", (WTraceReal)walker.Multiplicity  );
  bsr.collect("logpsi"      , (WTraceReal)wfn.getLogPsi()      );
  bsr.collect("phase"       , (WTraceReal)wfn.getPhase()       );
  //    from PropertyList
  if (bsr.first_collect)
  {
    for(size_t n=0;n<pset.PropertyList.size();++n)
    {
      auto& name  = pset.PropertyList.Names[n];
      auto& value = walker.Properties(0,n);
      if(properties_include.find(name) != properties_include.end())
      {
        bsr.collect(name, (WTraceReal)value );
        property_indices.push_back(n);
      }
      if(name=="LocalEnergy")
        energy_index = n;
    }
    if(energy_index<0)
      throw std::runtime_error("TraceCollector::collect  energy_index must not be negative");
  }
  else
    for(auto n: property_indices)
    {
      auto& name  = pset.PropertyList.Names[n];
      auto& value = walker.Properties(0,n);
      bsr.collect(name, (WTraceReal)value );
    }
  //    nodal proximity measures
  auto& gv = Gtmp.storage();
  WTraceReal dr = std::numeric_limits<WTraceReal>::max();
  for(size_t n=0; n<gv.size(); ++n)
    dr = std::min(dr,std::abs(1./std::real(gv[n])));
  auto dr_node_min = dr;
  WTraceReal dlogpsi2 = 0.;
  for(size_t n=0; n<gv.size(); ++n)
    dlogpsi2 += std::real(gv[n])*std::real(gv[n]);
  WTraceReal dphase2 = 0.;
  for(size_t n=0; n<gv.size(); ++n)
    dphase2 += std::imag(gv[n])*std::imag(gv[n]);
  bsr.collect("dr_node_min", dr_node_min);
  bsr.collect("dlogpsi2"   , dlogpsi2   ); // dr_node = 1/sqrt(dlogpsi2)
  bsr.collect("dphase2"    , dphase2    ); // dr_node = 1/sqrt(dphase2)  
  bsr.resetCollect();

  // save the energy of this walker
  steps.push_back((size_t)current_step);
  energies.push_back((WTraceReal)walker.Properties(0,energy_index));

}


void WalkerTraceCollector::resetBuffers()
{
  if (state.verbose) app_log() << "WalkerTraceCollector::reset_buffers"<<std::endl;
  walker_property_int_buffer.resetBuffer();
  walker_property_real_buffer.resetBuffer();
  walker_particle_real_buffer.resetBuffer();
  steps.resize(0);
  energies.resize(0);
}


void WalkerTraceCollector::checkBuffers()
{
  if (state.verbose) app_log() << "WalkerTraceCollector::checkBuffers"<<std::endl;
  size_t nrows = walker_property_int_buffer.nrows();
  auto prop_real_bad = walker_property_real_buffer.nrows()!=nrows;
  auto part_real_bad = walker_particle_real_buffer.nrows()!=nrows;
  auto steps_bad     = steps.size()!=nrows;
  auto energies_bad  = energies.size()!=nrows;
  auto any_bad = prop_real_bad || part_real_bad || steps_bad || energies_bad;

  if(prop_real_bad) app_log()<<"WalkerTraceCollector::checkBuffers  walker_property_real_buffer row count does not match\n";
  if(part_real_bad) app_log()<<"WalkerTraceCollector::checkBuffers  walker_particle_real_buffer row count does not match\n";
  if(steps_bad) app_log()<<"WalkerTraceCollector::checkBuffers  steps entry count does not match\n";
  if(energies_bad) app_log()<<"WalkerTraceCollector::checkBuffers  energies entry count does not match\n";
  if(any_bad)
    throw std::runtime_error("WalkerTraceCollector::checkBuffers  buffer lengths do not match");
}




////////////////////////////////
// WalkerTraceManager methods //
////////////////////////////////

WalkerTraceManager::WalkerTraceManager(WalkerTraceInput& inp, bool allow_traces, std::string series_root, Communicate* comm)
{
  state.reset();
  communicator              = comm;
  file_root                 = series_root;
  bool method_allows_traces = allow_traces;

  bool traces_requested     = inp.present;
  state.traces_active       = traces_requested && method_allows_traces;

  if (state.traces_active)
  {
    if (omp_get_thread_num() == 0)
    {
      app_log() << "\n  WalkerTraceManager::put() " << std::endl;
      app_log() << "    traces requested      : " << traces_requested << std::endl;
      app_log() << "    method allows traces  : " << method_allows_traces << std::endl;
      app_log() << "    traces active         : " << state.traces_active << std::endl;
      app_log() << std::endl;
    }
    state.step_period   = inp.get<int>("step_period");
    state.verbose       = inp.get<bool>("verbose");
    bool qtiles         = inp.get<bool>("qtiles");
    write_particle_data = inp.get<bool>("particle");
    write_min_data      = inp.get<bool>("min")    && qtiles;
    write_max_data      = inp.get<bool>("max")    && qtiles;
    write_med_data      = inp.get<bool>("median") && qtiles;
  }
  
  wmin_property_int_buffer.label  = "wmin_property_int";
  wmin_property_real_buffer.label = "wmin_property_real";
  wmin_particle_real_buffer.label = "wmin_particle_real";
  
  wmax_property_int_buffer.label  = "wmax_property_int";
  wmax_property_real_buffer.label = "wmax_property_real";
  wmax_particle_real_buffer.label = "wmax_particle_real";
  
  wmed_property_int_buffer.label  = "wmed_property_int";
  wmed_property_real_buffer.label = "wmed_property_real";
  wmed_particle_real_buffer.label = "wmed_particle_real";

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
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::startRun " << std::endl;
  for (auto& tc: collectors)
    tc->state = state;
  checkCollectors(collectors);
  openFile(collectors);
}


void WalkerTraceManager::stopRun()
{
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::stopRun " << std::endl;
  closeFile();
}


void WalkerTraceManager::writeBuffers(std::vector<WalkerTraceCollector*>& collectors)
{
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::writeBuffers "<<std::endl;

  if(write_min_data)
  {
    wmin_property_int_buffer.resetBuffer();
    wmin_property_real_buffer.resetBuffer();
    wmin_particle_real_buffer.resetBuffer();
  }
  if(write_max_data)
  {
    wmax_property_int_buffer.resetBuffer();
    wmax_property_real_buffer.resetBuffer();
    wmax_particle_real_buffer.resetBuffer();
  }
  if(write_med_data)
  {
    wmed_property_int_buffer.resetBuffer();
    wmed_property_real_buffer.resetBuffer();
    wmed_particle_real_buffer.resetBuffer();
  }

  // collect energy information and extract info from min/max/median energy walkers
  if(write_min_data || write_max_data || write_med_data)
  {
    for (size_t c=0; c<collectors.size(); ++c)
    {
      auto& tc = *collectors[c];
      tc.checkBuffers();
      for (size_t r=0; r<tc.energies.size(); ++r)
        energy_order.push_back(std::make_tuple(tc.steps[r],tc.energies[r],c,r));
    }
    std::sort(energy_order.begin(), energy_order.end());
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
  int nprocs = communicator->size();
  int rank   = communicator->rank();
  std::array<char, 32> ptoken;
  std::string file_name = file_root;
  if (nprocs > 1)
  {
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
  hdf_file        = std::make_unique<hdf_archive>();
  bool successful = hdf_file->create(file_name);
  if (!successful)
    throw std::runtime_error("WalkerTraceManager::openHDFFile  failed to open hdf file " + file_name);
}


void WalkerTraceManager::writeBuffersHDF(std::vector<WalkerTraceCollector*>& collectors)
{
  if (state.verbose) app_log() << "WalkerTraceManager::writeBuffersHDF " << std::endl;
  WalkerTraceCollector& tc_lead = *collectors[0];
  if(!registered_hdf)
  {
    tc_lead.walker_property_int_buffer.registerHDFData(*hdf_file);
    tc_lead.walker_property_real_buffer.registerHDFData(*hdf_file);
    if(write_particle_data)
      tc_lead.walker_particle_real_buffer.registerHDFData(*hdf_file);
    if(write_min_data)
    {
      wmin_property_int_buffer.registerHDFData(*hdf_file);
      wmin_property_real_buffer.registerHDFData(*hdf_file);
      wmin_particle_real_buffer.registerHDFData(*hdf_file);
    }
    if(write_max_data)
    {
      wmax_property_int_buffer.registerHDFData(*hdf_file);
      wmax_property_real_buffer.registerHDFData(*hdf_file);
      wmax_particle_real_buffer.registerHDFData(*hdf_file);
    }
    if(write_med_data)
    {
      wmed_property_int_buffer.registerHDFData(*hdf_file);
      wmed_property_real_buffer.registerHDFData(*hdf_file);
      wmed_particle_real_buffer.registerHDFData(*hdf_file);
    }
    registered_hdf = true;
  }
  for (int ip = 0; ip < collectors.size(); ++ip)
  {
    WalkerTraceCollector& tc = *collectors[ip];
    tc.walker_property_int_buffer.writeHDF(*hdf_file, tc_lead.walker_property_int_buffer.hdf_file_pointer);
    tc.walker_property_real_buffer.writeHDF(*hdf_file, tc_lead.walker_property_real_buffer.hdf_file_pointer);
    if(write_particle_data)
      tc.walker_particle_real_buffer.writeHDF(*hdf_file, tc_lead.walker_particle_real_buffer.hdf_file_pointer);
  }
  if(write_min_data)
  {
    wmin_property_int_buffer.writeHDF(*hdf_file);
    wmin_property_real_buffer.writeHDF(*hdf_file);
    wmin_particle_real_buffer.writeHDF(*hdf_file);
  }
  if(write_max_data)
  {
    wmax_property_int_buffer.writeHDF(*hdf_file);
    wmax_property_real_buffer.writeHDF(*hdf_file);
    wmax_particle_real_buffer.writeHDF(*hdf_file);
  }
  if(write_med_data)
  {
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
