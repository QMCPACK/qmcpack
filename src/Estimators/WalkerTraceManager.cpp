//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "WalkerTraceManager.h"

#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsData/AttributeSet.h"

#include "Particle/Walker.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"



namespace qmcplusplus
{

using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;


  
template<typename T>
void report_row(T& b)
{
  auto& buf = b.buffer;
  auto irow = buf.size(0)-1;
  for(size_t n=0;n<b.quantity_info.size();++n)
  {
    auto& q = b.quantity_info[n];
    auto  v = buf(irow,q.buffer_start);
    app_log()<<"  "<<q.name<<" = "<<v<<std::endl;
  }
}

template<typename T>
void report_buffer(WalkerTraceBuffer<T>& wtb)
{
  auto& b = wtb.buffer;
  for (size_t r=0;r<b.size(0);++r)
  {
    for (size_t c=0;c<b.size(1);++c)
      app_log()<<"  "<<b(r,c);
    app_log()<<"\n";
  }
}



///////////////////////////////////
// WalkerTraceCollector methods  //
///////////////////////////////////

WalkerTraceCollector::WalkerTraceCollector()
  : properties_include{"R2Accepted","R2Proposed","LocalEnergy","LocalPotential","Kinetic","ElecElec","ElecIon","LocalECP","NonLocalECP"}
{
  state.reset_permissions();
  energy_index = -1;
  energies.resize(0);
  walker_property_int_buffer.label  = "walker_property_int";
  walker_property_real_buffer.label = "walker_property_real";
  walker_particle_real_buffer.label = "walker_particle_real";
}



void WalkerTraceCollector::set_state(const WalkerTraceState& tms)
{
  state = tms;
  walker_property_int_buffer.verbose  = state.verbose;
  walker_property_real_buffer.verbose = state.verbose;
  walker_particle_real_buffer.verbose = state.verbose;
}


void WalkerTraceCollector::startBlock(int nsteps)
{
  if(!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceCollector::startBlock " << std::endl;
  reset_buffers();
}


void WalkerTraceCollector::collect(MCPWalker& walker, ParticleSet& pset, TrialWaveFunction& wfn, QMCHamiltonian& ham)
{
  if(!state.traces_active || pset.current_step%state.step_period!=0) return;

  //app_log()<<"TraceCollector::collect (step "<<pset.current_step<<")"<<std::endl;

  auto& bsi = walker_property_int_buffer;
  auto& bsr = walker_property_real_buffer;
  auto& bar = walker_particle_real_buffer;

  // collect integer walker properties
  bsi.collect("step"        , (WTraceInt)pset.current_step    );
  bsi.collect("id"          , (WTraceInt)walker.getWalkerID() );
  bsi.collect("parent_id"   , (WTraceInt)walker.getParentID() );
  bsi.collect("age"         , (WTraceInt)walker.Age           );
  bsi.reset_collect();

  // collect real walker properties
  bsr.collect("weight"      , (WTraceReal)walker.Weight        );
  bsr.collect("multiplicity", (WTraceReal)walker.Multiplicity  );
  bsr.collect("logpsi"      , (WTraceReal)wfn.getLogPsi()      );
  bsr.collect("phase"       , (WTraceReal)wfn.getPhase()       );
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
  bsr.reset_collect();
  
  // collect per particle walker quantities
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
  bar.reset_collect();

  // save the energy of this walker
  energies.push_back((WTraceReal)walker.Properties(0,energy_index));

  //app_log()<<"  walker_property_int buffer contents:\n";
  //report_buffer(bsi);

  //app_log()<<"\ninteger walker buffer:\n";
  //report_row(bsi);
  //app_log()<<"\nreal walker buffer:\n";
  //report_row(bsr);
  //
  //bsi.write_summary();
  //bsr.write_summary();
  //bar.write_summary();
  ////APP_ABORT("JTK");

}


void WalkerTraceCollector::reset_step() 
{ 
  energies.resize(0);
}


void WalkerTraceCollector::reset_buffers()
{
  if (state.verbose) app_log() << "WalkerTraceCollector::reset_buffers"<<std::endl;
  walker_property_int_buffer.reset_buffer();
  walker_property_real_buffer.reset_buffer();
  walker_particle_real_buffer.reset_buffer();
}




/////////////////////////////////
// WalkerTraceManager methods  //
/////////////////////////////////

WalkerTraceManager::WalkerTraceManager(Communicate* comm)
{
  state.reset_permissions();
  communicator   = comm;
  
  registered_hdf = false;
  
  // new walker buffers, etc
  write_particle_data = false;
  write_min_data      = true;
  write_max_data      = true;
  write_med_data      = true;
  
  wmin_property_int_buffer.label  = "wmin_property_int";
  wmin_property_real_buffer.label = "wmin_property_real";
  wmin_particle_real_buffer.label = "wmin_particle_real";
  
  wmax_property_int_buffer.label  = "wmax_property_int";
  wmax_property_real_buffer.label = "wmax_property_real";
  wmax_particle_real_buffer.label = "wmax_particle_real";
  
  wmed_property_int_buffer.label  = "wmed_property_int";
  wmed_property_real_buffer.label = "wmed_property_real";
  wmed_particle_real_buffer.label = "wmed_particle_real";
}


void WalkerTraceManager::put(xmlNodePtr cur, bool allow_traces, std::string series_root)
{
  if (state.verbose) app_log()<<"WalkerTraceManager::put"<<std::endl;
  state.reset_permissions();
  bool traces_requested     = cur != NULL;
  bool method_allows_traces = allow_traces;
  file_root                 = series_root;
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
    //read trace attributes
    std::string verbose_write = "no";
    OhmmsAttributeSet attrib;
    attrib.add(state.step_period, "step_period");
    attrib.add(verbose_write, "verbose");
    attrib.put(cur);
    state.verbose = verbose_write == "yes";
  }
}


WalkerTraceCollector* WalkerTraceManager::makeCollector()
{
  if (state.verbose) app_log() << "WalkerTraceManager::makeCollector " << std::endl;
  WalkerTraceCollector* tc = new WalkerTraceCollector();
  tc->set_state(get_state());
  return tc;
}


void WalkerTraceManager::startRun(int blocks, std::vector<WalkerTraceCollector*>& collectors)
{
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::startRun " << std::endl;
  check_collectors(collectors);
  open_file(collectors);
}


void WalkerTraceManager::stopRun()
{
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::stopRun " << std::endl;
  close_file();
}


void WalkerTraceManager::stopStep()
{
  if(!state.traces_active) return;
  // find min/max/median walker and collect in buffers
}


void WalkerTraceManager::write_buffers(std::vector<WalkerTraceCollector*>& collectors, int block)
{
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::write_buffers "<<std::endl;
  write_buffers_hdf(collectors);
}


void WalkerTraceManager::check_collectors(std::vector<WalkerTraceCollector*>& collectors)
{
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::check_collectors" << std::endl;
  if (collectors.size() > 0)
  {
    bool all_same = true;
    WalkerTraceCollector& ref = *collectors[0];
    for (int i = 0; i < collectors.size(); ++i)
    {
      WalkerTraceCollector& tm = *collectors[i];
      all_same &= tm.walker_property_int_buffer.same_as(ref.walker_property_int_buffer);
      all_same &= tm.walker_property_real_buffer.same_as(ref.walker_property_real_buffer);
      all_same &= tm.walker_particle_real_buffer.same_as(ref.walker_particle_real_buffer);
    }
    if (!all_same)
    {
      throw std::runtime_error("WalkerTraceManager::check_collectors  trace buffer widths of collectors do not match\n  contiguous write is "
                               "impossible\n  this was first caused by collectors contributing array traces from identical, but "
                               "differently named, particlesets such as e, e2, e3 ... (fixed)\n  please check the WalkerTraceManager "
                               "summaries printed above");
    }
  }
}


void WalkerTraceManager::open_file(std::vector<WalkerTraceCollector*>& collectors)
{
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::open_file "<<std::endl;
  open_hdf_file(collectors);
}


void WalkerTraceManager::close_file()
{
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::close_file " << std::endl;
  close_hdf_file();
}


void WalkerTraceManager::open_hdf_file(std::vector<WalkerTraceCollector*>& collectors)
{
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::open_hdf_file " << std::endl;
  if (collectors.size() == 0) 
    throw std::runtime_error("WalkerTraceManager::open_hdf_file  no trace collectors exist, cannot open file");
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
  if (state.verbose) app_log() << "WalkerTraceManager::open_hdf_file  opening traces hdf file " << file_name << std::endl;
  hdf_file        = std::make_unique<hdf_archive>();
  bool successful = hdf_file->create(file_name);
  if (!successful)
    throw std::runtime_error("WalkerTraceManager::open_hdf_file  failed to open hdf file " + file_name);
}


void WalkerTraceManager::write_buffers_hdf(std::vector<WalkerTraceCollector*>& collectors)
{
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::write_buffers_hdf " << std::endl;
  WalkerTraceCollector& tc_lead = *collectors[0];
  if(!registered_hdf)
  {
    tc_lead.walker_property_int_buffer.register_hdf_data(*hdf_file);
    tc_lead.walker_property_real_buffer.register_hdf_data(*hdf_file);
    tc_lead.walker_particle_real_buffer.register_hdf_data(*hdf_file);
    registered_hdf = true;
  }
  for (int ip = 0; ip < collectors.size(); ++ip)
  {
    WalkerTraceCollector& tc = *collectors[ip];
    tc.walker_property_int_buffer.write_hdf(*hdf_file, tc_lead.walker_property_int_buffer.hdf_file_pointer);
    tc.walker_property_real_buffer.write_hdf(*hdf_file, tc_lead.walker_property_real_buffer.hdf_file_pointer);
    tc.walker_particle_real_buffer.write_hdf(*hdf_file, tc_lead.walker_particle_real_buffer.hdf_file_pointer);
  }
}


void WalkerTraceManager::close_hdf_file()
{ 
  if (!state.traces_active) return;
  if (state.verbose) app_log() << "WalkerTraceManager::close_hdf_file " << std::endl;
  hdf_file.reset(); 
}



} // namespace qmcplusplus
