#ifndef QMCPLUSPLUS_AFQMC_ENERGYESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_ENERGYESTIMATOR_H

#include<Message/MPIObjectBase.h>
#include"AFQMC/config.h"
#include<vector>
#include<queue>
#include<string>
#include<iostream>
#include<fstream>

#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/NewTimer.h"

#include "boost/multi_array.hpp"

#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"

namespace qmcplusplus
{
namespace afqmc
{

class EnergyEstimator: public EstimatorBase 
{

  public:

  EnergyEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo info, xmlNodePtr cur, 
        Wavefunction& wfn, bool timer=true):
            EstimatorBase(info),TG(tg_),wfn0(wfn)
  {

    data.resize(2);
  }

  ~EnergyEstimator() {}

  void accumulate_step(WalkerSet& wlks, std::vector<ComplexType>& curData) {}

  void accumulate_block(WalkerSet& wset)
  {
    AFQMCTimers[energy_timer]->start();
    size_t nwalk = wset.size();
    if(eloc.shape()[0] != nwalk || eloc.shape()[1] != 3)
      eloc.resize(boost::extents[nwalk][3]);
    if(ovlp.shape()[0] != nwalk) 
      ovlp.resize(boost::extents[nwalk]);

    ComplexType dum, et;
    wfn0.Energy(wset,eloc,ovlp);
    if(TG.TG_local().root()) {
      data[0] = data[1] = std::complex<double>(0,0);  
      for(int i=0; i<nwalk; i++) {
        auto wi = wset[i];
        if(std::isnan(real(wi.weight()))) continue;
        dum = wi.weight()*ovlp[i]/wi.overlap();
        et = eloc[i][0]+eloc[i][1]+eloc[i][2];
        if( (!std::isfinite(real(dum))) || (!std::isfinite(real(et*dum))) ) continue; 
        data[1] += dum; 
        data[0] += et*dum;
      }
      TG.TG_heads().all_reduce_in_place_n(data.begin(),data.size(),std::plus<>());
    } 
    AFQMCTimers[energy_timer]->stop();

  }

  void tags(std::ofstream& out) 
  {
    if(TG.Global().root()) {
      out<<"EnergyEstim_" <<name <<"_nume_real  EnergyEstim_" <<name <<"_nume_imag " 
         <<"EnergyEstim_" <<name <<"_deno_real  EnergyEstim_" <<name <<"_deno_imag " 
         <<"EnergyEstim_" <<name <<"_timer ";
    }
  }

  void print(std::ofstream& out,WalkerSet& wset)
  {
    if(TG.Global().root()) {
     int n = wset.get_global_target_population(); 
      out<< data[0].real()/n << " " << data[0].imag()/n << " " 
         << data[1].real()/n << " " << data[1].imag()/n << " "
         <<AFQMCTimers[energy_timer]->get_total() <<" "; 
      AFQMCTimers[energy_timer]->reset();
    }
  }

  private:

  std::string name;

  TaskGroup_& TG;

  Wavefunction& wfn0;

  boost::multi_array<ComplexType,2> eloc;
  boost::multi_array<ComplexType,1> ovlp;

  std::vector<std::complex<double> > data;

};
}
}

#endif
