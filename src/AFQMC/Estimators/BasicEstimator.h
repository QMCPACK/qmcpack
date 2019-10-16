#ifndef QMCPLUSPLUS_AFQMC_BASICESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_BASICESTIMATOR_H

#include"AFQMC/config.h"
#include<vector>
#include<queue>
#include<string>
#include<iostream>
#include<fstream>

#include "OhmmsData/libxmldefs.h"

#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"


namespace qmcplusplus
{
namespace afqmc
{

class BasicEstimator: public EstimatorBase
{

  public:

  BasicEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo info,
        std::string title, xmlNodePtr cur, bool impsamp_):
                EstimatorBase(info),TG(tg_),
                nwfacts(0),writer(false),importanceSampling(impsamp_),
                timers(false)
  {

    if(cur!=NULL)  {
      ParameterSet m_param;
      std::string str2;
      m_param.add(str2,"timers","std::string");
      m_param.add(nwfacts,"nhist","int");
      m_param.put(cur);

      std::transform(str2.begin(),str2.end(),str2.begin(),(int (*)(int)) tolower);
      if(str2 == "yes" || str2 == "true")
        timers = true;
    }

#ifndef ENABLE_TIMERS
    timers=false;
#endif

    assert(nwfacts>=0);
    weight_product=ComplexType(1.0,0.0);
    for(int i=0; i<nwfacts; i++)
      weight_factors.push(weight_product);

    app_log()<<"  BasicEstimator: Number of products in weight history: " <<nwfacts <<std::endl;

    writer = (TG.getGlobalRank()==0);

    data.resize(10);
    data2.resize(10);
    data3.resize(2);

    if(timers) {
      AFQMCTimers[block_timer]->reset();
      AFQMCTimers[pseudo_energy_timer]->reset();
      AFQMCTimers[vHS_timer]->reset();
      AFQMCTimers[vbias_timer]->reset();
      AFQMCTimers[G_for_vbias_timer]->reset();
      AFQMCTimers[propagate_timer]->reset();
      AFQMCTimers[E_comm_overhead_timer]->reset();
      AFQMCTimers[vHS_comm_overhead_timer]->reset();
      AFQMCTimers[assemble_X_timer]->reset();
      AFQMCTimers[popcont_timer]->reset();
      AFQMCTimers[ortho_timer]->reset();
      AFQMCTimers[setup_timer]->reset();
      AFQMCTimers[extra_timer]->reset();
    }

    enume=0.0;
    edeno=0.0;
    enume_sub=0.0;
    edeno_sub=0.0;
    enume2=0.0;
    edeno2=0.0;
    weight=0.0;
    weight_sub=0.0;
    nwalk = 0;
    nwalk_good = 0;
    nwalk_sub = 0;
    ncalls=0;
    ncalls_substep=0;
    nwalk_min=1000000;
    nwalk_max=0;

    // first block will always be off, move to Driver if problematic
    AFQMCTimers[block_timer]->start();

  }

  ~BasicEstimator() {}

  void accumulate_block(WalkerSet& wset) {}


  //  curData:
  //  0: inverse of the factor used to rescale the weights
  //  1: 1/nW * sum_i w_i * Eloc_i   (where w_i is the normalized weight)
  //  2: 1/nW * sum_i w_i            (where w_i is the normalized weight)
  //  3: sum_i abs(w_i)       (where w_i is the normalized weight)
  //  4: 1/nW * sum_i abs(<psi_T|phi_i>)
  //  5: nW                          (total number of walkers)
  //  6: "healthy" nW                (total number of "healthy" walkers)
  void accumulate_step(WalkerSet& wset, std::vector<ComplexType>& curData)
  {

    ncalls++;
    if(nwfacts>0) {
      weight_product *= (curData[0]/weight_factors.front());
      weight_factors.pop();
      weight_factors.push(curData[0]);
    } else
      weight_product = ComplexType(1.0,0.0);

    data2[0] = curData[1].real();
    data2[1] = curData[2].real();

    int nwlk = wset.size();
    if(nwlk>nwalk_max) nwalk_max=nwlk;
    if(nwlk<nwalk_min) nwalk_min=nwlk;
    enume += (curData[1]/curData[2])*weight_product;
    edeno += weight_product;
    weight += curData[3].real();
    ovlp += wset.getLogOverlapFactor(); //curData[4].real();
    nwalk += static_cast<int>(std::floor(curData[5].real()));
    nwalk_good += static_cast<int>(std::floor(curData[6].real()));
  }

  void tags(std::ofstream& out)
  {
    if(writer) {
      if(nwfacts>0) {
        out<<"nWalkers weight Eloc_nume Eloc_deno ";
      } else {
        out<<"nWalkers weight PseudoEloc ";
      }
      out<<"LogOvlpFactor ";
    }
  }

  void tags_timers(std::ofstream& out)
  {
    if(writer)
      if(timers) out<<"PseudoEnergy_t vHS_t vbias_t G_t Propagate_t Energy_comm_t vHS_comm_t X_t popC_t ortho_t setup_t extra_t Block_t ";
  }

  void print(std::ofstream& out, hdf_archive& dump, WalkerSet& wset)
  {
    data[0] = enume.real()/ncalls;
    data[1] = edeno.real()/ncalls;

    if(writer) {
      out<<std::setprecision(6) <<nwalk/ncalls <<" " <<weight/ncalls <<" " <<std::setprecision(16);
      if(nwfacts>0) {
        out<<enume.real()/ncalls <<" " <<edeno.real()/ncalls <<" ";
      } else {
        out<<enume.real()/ncalls <<" ";
      }
      out<<ovlp/ncalls <<" ";
    }

    enume=0.0;
    edeno=0.0;
    weight=0.0;
    enume2=0.0;
    edeno2=0.0;
    ncalls=0;
    nwalk=0;
    nwalk_good=0;
    nwalk_min=1000000;
    nwalk_max=0;
    ovlp=0;
  }

  void print_timers(std::ofstream& out)
  {
    AFQMCTimers[block_timer]->stop();

    if(writer) {
      if(timers) out<<std::setprecision(5)
                    <<AFQMCTimers[pseudo_energy_timer]->get_total() <<" "
                    <<AFQMCTimers[vHS_timer]->get_total() <<" "
                    <<AFQMCTimers[vbias_timer]->get_total() <<" "
                    <<AFQMCTimers[G_for_vbias_timer]->get_total() <<" "
                    <<AFQMCTimers[propagate_timer]->get_total() <<" "
                    <<AFQMCTimers[E_comm_overhead_timer]->get_total() <<" "
                    <<AFQMCTimers[vHS_comm_overhead_timer]->get_total() <<" "
                    <<AFQMCTimers[assemble_X_timer]->get_total() <<" "
                    <<AFQMCTimers[popcont_timer]->get_total() <<" "
                    <<AFQMCTimers[ortho_timer]->get_total() <<" "
                    <<AFQMCTimers[setup_timer]->get_total() <<" "
                    <<AFQMCTimers[extra_timer]->get_total() <<" "
                    <<AFQMCTimers[block_timer]->get_total() <<" "
                    <<std::setprecision(16);
    }
    if(timers) {
      AFQMCTimers[block_timer]->reset();
      AFQMCTimers[pseudo_energy_timer]->reset();
      AFQMCTimers[vHS_timer]->reset();
      AFQMCTimers[vbias_timer]->reset();
      AFQMCTimers[G_for_vbias_timer]->reset();
      AFQMCTimers[propagate_timer]->reset();
      AFQMCTimers[E_comm_overhead_timer]->reset();
      AFQMCTimers[vHS_comm_overhead_timer]->reset();
      AFQMCTimers[popcont_timer]->reset();
      AFQMCTimers[assemble_X_timer]->reset();
      AFQMCTimers[ortho_timer]->reset();
      AFQMCTimers[setup_timer]->reset();
      AFQMCTimers[extra_timer]->reset();
    }
    AFQMCTimers[block_timer]->start();
  }

  double getEloc()
  {
    return data[0]/data[1];
  }

  double getEloc_step()
  {
    return data2[0]/data2[1];
  }


  private:

  afqmc::TaskGroup_& TG;

  int nwfacts;

  bool writer;

  bool importanceSampling;

  std::vector<double> data, data2, data3;

  std::queue<ComplexType> weight_factors;
  ComplexType weight_product=ComplexType(1.0,0.0);

  ComplexType enume=0.0,edeno=0.0;
  ComplexType enume_sub=0.0,edeno_sub=0.0;
  ComplexType enume2=0.0,edeno2=0.0;
  RealType weight, weight_sub, ovlp, ovlp_sub;
  int nwalk_good, nwalk, ncalls, ncalls_substep, nwalk_sub, nwalk_min, nwalk_max;

  // optional
  bool timers;

};
}
}

#endif
