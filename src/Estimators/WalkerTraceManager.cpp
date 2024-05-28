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

#include "Particle/Walker.h"
//#include "Particle/ParticleSet.h"
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
    {
      app_log()<<"  "<<b(r,c);
    }
    app_log()<<"\n";
  }
}




void WalkerTraceCollector::collect(MCPWalker& walker, ParticleSet& pset, TrialWaveFunction& wfn, QMCHamiltonian& ham)
{
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

}
