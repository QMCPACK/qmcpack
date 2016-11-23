#ifndef QMCPLUSPLUS_AFQMC_BASICESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_BASICESTIMATOR_H

#include<Message/MPIObjectBase.h>
#include"AFQMC/config.h"
#include<vector>
#include<string>
#include<iostream>
#include<fstream>

#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Wavefunctions/WavefunctionHandler.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"
#include "AFQMC/Propagators/PropagatorBase.h"
#include "AFQMC/Estimators/SlaterDetOperations.h"

namespace qmcplusplus
{

class BasicEstimator: public EstimatorBase 
{

  public:

  typedef HamiltonianBase* HamPtr;
  typedef WavefunctionHandler* WfnPtr;
  typedef PropagatorBase* PropPtr;
  typedef WalkerHandlerBase* WSetPtr;

  BasicEstimator(Communicate *c):EstimatorBase(c),EstimEloc(false), timers(true), prtnwalk(true),
                                 overlap(false), diag(false), nstates(0),SDet(c) 
  {}

  ~BasicEstimator() {}


  void accumulate_block(WSetPtr wlkBucket)
  {

    int nW = wlkBucket->numWalkers();
    ComplexType *sm,eloc,dum,oa,ob,w,ooa,oob;
    LocalTimer->start("Block::EstimatorEloc");

    if(EstimEloc) {
      wfn0->evaluateLocalEnergyAndOverlap("Estimator",-1,wlkBucket);

      if(core_rank==0) {
        for(int i=0, cnt=0; i<nW; i++) {
          if(!wlkBucket->isAlive(i) || std::abs(wlkBucket->getWeight(i)) <= 1e-6 || std::isnan(wlkBucket->getWeight(i).real())) continue;
      
          sm = wlkBucket->getWalker(i,w,dum,ooa,oob);
          if(std::isnan(ooa.real()) || std::isnan(oob.real()) || std::abs(oa*ob) < 1e-8  || std::abs(ooa*oob) < 1e-8) w=0;
          dum = w*oa*ob/(ooa*oob);
          edeno2 += dum; 
          enume2 += eloc*dum;
        }
      }
    }
    LocalTimer->stop("Block::EstimatorEloc");
  }


  void accumulate_step(WSetPtr wlkBucket)
  {

    ncalls++;
    int nw=0;
    RealType instant_weight=0.0;
    if(core_rank==0) {
      int nW = wlkBucket->numWalkers();
      enume_sub = edeno_sub = 0.0;
      ComplexType w,oa,ob,eloc;
      for(int i=0; i<nW; i++) {
        ComplexType* dum = wlkBucket->getWalker(i,w,eloc,oa,ob);
        if(!wlkBucket->isAlive(i) || std::abs(w) <= 1e-6) continue;
        enume_sub += w*eloc;
        edeno_sub += w;
        nw++;
        instant_weight += std::abs(w);
      }
      if(nw>nwalk_max) nwalk_max=nw;
      if(nw<nwalk_min) nwalk_min=nw;

      data2[0] = enume_sub.real();
      data2[1] = edeno_sub.real();
      data2[4] = nwalk_sub/ncalls_substep;
      data2[5] = instant_weight; 
      data2[6] = weight_sub/ncalls_substep;
      data2[7] = ovlp_sub/ncalls_substep;

      myComm->allreduce(data2,MPI_COMM_TG_LOCAL_HEADS);

      wlkBucket->scaleWeight(targetW/data2[5]);

      data2[0] *= targetW/data2[5];
      data2[1] *= targetW/data2[5];

      nwalk += nwalk_sub/ncalls_substep;
      enume += enume_sub*targetW/data2[5];
      edeno += edeno_sub*targetW/data2[5];
      weight += instant_weight; 
      ovlp += ovlp_sub/ncalls_substep;

      enume_sub=0.0;
      edeno_sub=0.0;
      nwalk_sub = 0;
      weight_sub=0.0;
      ncalls_substep=0;
      ovlp_sub=0;

    }   

    myComm->bcast(data2,MPI_COMM_TG_LOCAL);

  }

  void accumulate_substep(WSetPtr wlkBucket)
  {

    if(core_rank != 0) return;
    ncalls_substep++;
    int nW = wlkBucket->numWalkers();
    ComplexType w,oa,ob,eloc;
    int cnt1=0; 
    RealType sumo=0;
    for(int i=0; i<nW; i++) {
      ComplexType* dum = wlkBucket->getWalker(i,w,eloc,oa,ob);
      if(!wlkBucket->isAlive(i) || std::abs(w) <= 1e-6) continue;
      cnt1++;
      nwalk_sub++;
      enume_sub += w*eloc;
      edeno_sub += w;
      weight_sub += std::abs(w);
      sumo += std::abs(oa*ob);
    }
    if(cnt1>1) ovlp_sub+=sumo/static_cast<double>(cnt1);

  }

  void tags(std::ofstream& out) 
  {
    if(myComm->rank() == 0) {
      out<<"nWalkers weight  Eloc_nume Eloc_deno  "; 
      if(EstimEloc) out<<"ElocEstim_nume ElocEstim_deno  ";
      out<<"Ovlp "; 
      if(diag && nstates>0) {
        for(int i=0; i<nstates; i++) out<<"Ediag_" <<i <<" ";
#ifdef AFQMC_TIMER
        out<<"evaluate_H_S solve_GEV" <<" ";
#endif
      }
      if(timers && EstimEloc) out<<"TimeEstimEloc  "; 
      if(timers) out<<"TimePropg  TimePopControl TimeLoadBalance TimeOrtho TimeCommSetup TimeCommResize TimeCommExch TimeBlock "; 
      if(prtnwalk) out<<"nWmin nWmax ";
    }
  }

  void print(std::ofstream& out,WSetPtr wlkBucket)
  {

    // average over processors
    if(core_rank==0) {
      data[0] = enume.real()/ncalls;
      data[1] = edeno.real()/ncalls;
      data[2] = enume2.real();
      data[3] = edeno2.real();
      data[4] = nwalk/ncalls;
      data[6] = weight/ncalls;
      data[7] = ovlp/num_heads_tg;

      data[5] = 0.0;
      int nW = wlkBucket->numWalkers();
      for(int i=0; i<nW; i++) 
       if(wlkBucket->isAlive(i)) data[5] += std::abs(wlkBucket->getWeight(i));
      myComm->allreduce(data,MPI_COMM_TG_LOCAL_HEADS);
    }

    myComm->bcast(data,MPI_COMM_TG_LOCAL);

    if(writer) {
      out<<std::setprecision(6) <<data[4] <<" " <<data[6] <<" " <<std::setprecision(16) <<data[0] <<" " <<data[1] <<" ";
      if(EstimEloc) out<<data[2] <<" " <<data[3] <<" ";
      out<<data[7] <<" ";
      if(diag && nstates>0) {
        for(int i=0; i<nstates; i++) out<<eigVal[i] <<" ";
#ifdef AFQMC_TIMER
        out<<LocalTimer->total("SlaterDetOperations::diag::evaluate_H_S") <<" "
           <<LocalTimer->total("SlaterDetOperations::diag::solve_GEV") <<" ";
        LocalTimer->reset("SlaterDetOperations::diag::evaluate_H_S");
        LocalTimer->reset("SlaterDetOperations::diag::solve_GEV");
#endif
      }
      if(timers && EstimEloc)
                 out<<std::setprecision(3) <<LocalTimer->total("Block::EstimatorEloc") <<" ";
      if(timers) out<<std::setprecision(3) <<LocalTimer->total("SubStep::Propagate") <<" "
                    <<LocalTimer->total("Step::PopControl") <<" " 
                    <<LocalTimer->total("Step::loadBalance") <<" " 
                    <<LocalTimer->total("Step::Orthogonalize") <<" " 
                    <<LocalTimer->total("WalkerHandler::loadBalance::setup") <<" " 
                    <<LocalTimer->total("WalkerHandler::loadBalance::resize") <<" " 
                    <<LocalTimer->total("WalkerHandler::loadBalance::exchange") <<" " 
                    <<LocalTimer->total("Block::TOTAL") <<" "; 
      if(prtnwalk) out<<wlkBucket->nwalk_min <<" " <<wlkBucket->nwalk_max <<std::setprecision(12) <<" ";
    }

    enume=0.0;
    edeno=0.0;
    weight=0.0; 
    enume2=0.0;
    edeno2=0.0;
    ncalls=0;
    nwalk=0;
    nwalk_min=1000000;
    nwalk_max=0;
    ovlp=0;

    LocalTimer->reset("SubStep::Propagate");
    LocalTimer->reset("Block::EstimatorEloc");
    LocalTimer->reset("Step::PopControl"); 
    LocalTimer->reset("Step::loadBalance"); 
    LocalTimer->reset("Step::Orthogonalize"); 
    LocalTimer->reset("WalkerHandler::loadBalance::setup"); 
    LocalTimer->reset("WalkerHandler::loadBalance::exchange"); 
    LocalTimer->reset("WalkerHandler::loadBalance::resize"); 
    LocalTimer->reset("Block::TOTAL"); 
    
  }

  double getWeight()
  {
    return data[5]; 
  }

  double getWeight_step()
  {
    return data2[5]; 
  }

  double getEloc()
  {
    return data[0]/data[1]; 
  }

  double getEloc_step()
  {
    return data2[0]/data2[1]; 
  }

  void average(WSetPtr wlks) {} 

  bool parse(xmlNodePtr cur) 
  {
    if(cur==NULL) return false;
  
    ParameterSet m_param;
    std::string str1,str2,str3,str4,str5;
    m_param.add(str1,"estim_eloc","std::string");
    m_param.add(str1,"estimeloc","std::string");
    m_param.add(str1,"EstimEloc","std::string");
    m_param.add(str2,"timers","std::string");
    m_param.add(str3,"nwalk","std::string");
    m_param.add(str4,"overlap","std::string");
    m_param.add(str5,"diag","std::string");
    m_param.add(nstates,"nstates","int");
    m_param.put(cur);

    std::transform(str1.begin(),str1.end(),str1.begin(),(int (*)(int)) tolower);
    std::transform(str2.begin(),str2.end(),str2.begin(),(int (*)(int)) tolower);
    std::transform(str3.begin(),str3.end(),str3.begin(),(int (*)(int)) tolower);
    std::transform(str4.begin(),str4.end(),str4.begin(),(int (*)(int)) tolower);
    std::transform(str5.begin(),str5.end(),str5.begin(),(int (*)(int)) tolower);
    if(str1 == "yes" || str1 == "true")
      EstimEloc = true;
    if(str2 == "no" || str2 == "false")
      timers = false;
    if(str3 == "no" || str3 == "false")
      prtnwalk = false;
    if(str1 != "") EstimEloc_present=true; 
    if(str4 == "yes" || str4 == "true") overlap=true;
    if(str5 == "yes" || str5 == "true") diag=true;
    
    return true;
  }

  bool setup(std::vector<int>& TGdata, ComplexSMVector *v, HamiltonianBase* ham, WavefunctionHandler* wfn,myTimer* timer, MPI_Comm heads_comm, MPI_Comm tg_comm, MPI_Comm node_comm, MPI_Comm heads_of_tg_comm) 
  {
    ncores_per_TG=TGdata[4];
    core_rank = TGdata[1]%ncores_per_TG; 
    ham0=ham;
    wfn0=wfn;
    if(!EstimEloc_present) EstimEloc = wfn0->check_initialized("Estimator");
    writer = (myComm->rank()==0);
    LocalTimer = timer;

    MPI_COMM_TG_LOCAL_HEADS = heads_of_tg_comm;
    MPI_COMM_HEAD_OF_NODES = heads_comm;
    MPI_COMM_NODE_LOCAL = node_comm;
    MPI_COMM_TG_LOCAL = tg_comm;
    MPI_Comm_size(MPI_COMM_TG_LOCAL_HEADS,&num_heads_tg);

    if(overlap || (diag && nstates>0) ) {
      SDet.copyInfo(*this);
      SDet.setup(ham,LocalTimer);
      eigVal.resize(nstates);
      eigVec.resize(nstates,nstates);
    }
 
    data.resize(10);
    data2.resize(10);

    enume=0.0;
    edeno=0.0;
    enume_sub=0.0;
    edeno_sub=0.0;
    enume2=0.0;
    edeno2=0.0;
    weight=0.0;
    weight_sub=0.0;
    nwalk = 0;
    nwalk_sub = 0;
    ncalls=0;
    ncalls_substep=0;
    nwalk_min=1000000;
    nwalk_max=0;

    return true;
  }
  
  void setTargetWeight(RealType w0) { targetW = w0; } 

  private:

  std::vector<double> data, data2;  

  ComplexType enume=0.0,edeno=0.0;
  ComplexType enume_sub=0.0,edeno_sub=0.0;
  ComplexType enume2=0.0,edeno2=0.0;
  RealType weight, weight_sub, ovlp, ovlp_sub;
  RealType targetW=1;
  int nwalk, ncalls, ncalls_substep, nwalk_sub, nwalk_min, nwalk_max; 
  int core_rank;
  int ncores_per_TG;

  // optional 
  bool EstimEloc, timers, prtnwalk;    
  bool EstimEloc_present;

  bool overlap, diag;
  int nstates; 
  SlaterDetOperations SDet;

  std::vector<RealType> eigVal;  
  ComplexMatrix eigVec;
  ComplexType exactEnergy;

  myTimer* LocalTimer;

  int num_heads_tg;
  MPI_Comm MPI_COMM_HEAD_OF_NODES;
  MPI_Comm MPI_COMM_NODE_LOCAL;
  MPI_Comm MPI_COMM_TG_LOCAL;
  MPI_Comm MPI_COMM_TG_LOCAL_HEADS;

};
}

#endif
