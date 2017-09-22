#ifndef QMCPLUSPLUS_AFQMC_BASICESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_BASICESTIMATOR_H

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

#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Wavefunctions/WavefunctionHandler.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"
#include "AFQMC/Estimators/SlaterDetOperations.h"

namespace qmcplusplus
{

class BasicEstimator: public EstimatorBase 
{

  public:

  typedef HamiltonianBase* HamPtr;
  typedef WavefunctionHandler* WfnPtr;
  typedef WalkerHandlerBase* WSetPtr;

  BasicEstimator(Communicate *c):EstimatorBase(c),EstimEloc(false), timers(true), prtnwalk(false),
                                 overlap(false), diag(false), nstates(0),SDet(c),EstimEloc_present(false) 
  {}

  ~BasicEstimator() {}

  void accumulate_block(WSetPtr wlkBucket)
  {

    int nW = wlkBucket->numWalkers(true);
    ComplexType *sm,eloc,dum,oa,ob,w,ooa,oob;
    LocalTimer->start("Block::EstimatorEloc");
    BlockTimer->start();

    if(EstimEloc) {
      wfn0->evaluateLocalEnergyAndOverlap("Estimator",-1,wlkBucket);
      if(core_rank==0) {
        for(int i=0, cnt=0; i<nW; i++) {
          if(!wlkBucket->isAlive(i) || std::abs(wlkBucket->getWeight(i)) <= 1e-6 || std::isnan(wlkBucket->getWeight(i).real())) continue;
          sm = wlkBucket->getWalker(i,w,dum,ooa,oob);  // "impsampl"
          sm = wlkBucket->getWalker2(i,eloc,oa,ob);     // "estimator"
          dum = weight_product*w*oa*ob/(ooa*oob);
          //if( (!std::isfinite(dum.real())) || (!std::isfinite((eloc*dum).real())) || std::abs(oa*ob) < 1e-8  || std::abs(ooa*oob) < 1e-8) continue; 
          if( (!std::isfinite(dum.real())) || (!std::isfinite((eloc*dum).real())) ) continue; 
          edeno2 += dum; 
          enume2 += eloc*dum;
        }
        data3[0] = enume2.real();
        data3[1] = edeno2.real();
        myComm->allreduce(data3,MPI_COMM_TG_LOCAL_HEADS);
        data3[0] *= weight_product.real()/data3[1]; 
        data3[1] = weight_product.real(); 
      }
    }
    LocalTimer->stop("Block::EstimatorEloc");
    BlockTimer->stop();
  }


  //  curData:
  //  0: inverse of the factor used to rescale the weights
  //  1: 1/nW * sum_i w_i * Eloc_i   (where w_i is the normalized weight)
  //  2: 1/nW * sum_i w_i            (where w_i is the normalized weight)
  //  3: sum_i abs(w_i)       (where w_i is the normalized weight)
  //  4: 1/nW * sum_i abs(<psi_T|phi_i>)
  //  5: nW                          (total number of walkers)  
  //  6: "healthy" nW                (total number of "healthy" walkers)  
  void accumulate_step(WSetPtr wlkBucket, std::vector<ComplexType>& curData)
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

    int nwlk = wlkBucket->numWalkers(false);
    if(nwlk>nwalk_max) nwalk_max=nwlk;
    if(nwlk<nwalk_min) nwalk_min=nwlk;
    enume += (curData[1]/curData[2])*weight_product; 
    edeno += weight_product;
    weight += curData[3].real(); 
    ovlp += curData[4].real();
    nwalk += static_cast<int>(std::floor(curData[5].real())); 
    nwalk_good += static_cast<int>(std::floor(curData[6].real())); 
  }

  void tags(std::ofstream& out) 
  {
    if(myComm->rank() == 0) {
      if(nwfacts>0) {
        out<<"nWalkers weight  Eloc_nume Eloc_deno  "; 
        if(EstimEloc) out<<"ElocEstim_nume ElocEstim_deno  ";
      } else {
        out<<"nWalkers weight  Eloc  "; 
        if(EstimEloc) out<<"ElocEstim ";
      }
      out<<"Ovlp "; 
      if(diag && nstates>0) {
        for(int i=0; i<nstates; i++) out<<"Ediag_" <<i <<" ";
#ifdef AFQMC_TIMER
        out<<"evaluate_H_S solve_GEV" <<" ";
#endif
      }
      if(timers && EstimEloc) out<<"TimeEstimEloc  "; 
      if(timers) out<<"TimePropg  TimePopControl TimeOrtho TimeBranching TimeIdle TimeCommExch MaxTimeCommExch TimeBlock "; 
      if(prtnwalk) out<<"MaxBranch MaxExch ";
    }
  }

  void print(std::ofstream& out,WSetPtr wlkBucket)
  {

    data[0] = enume.real()/ncalls;
    data[1] = edeno.real()/ncalls;    

    double max_exch_time=0;
    if(timers) { 
      double t = LocalTimer->total("WalkerHandler::loadBalance::exchange");
      MPI_Reduce(&t,&max_exch_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    }

    if(writer) {
      out<<std::setprecision(6) <<nwalk/ncalls <<" " <<weight/ncalls <<" " <<std::setprecision(16);
      if(nwfacts>0) {
        out<<enume.real()/ncalls <<" " <<edeno.real()/ncalls <<" ";
        if(EstimEloc) out<<data3[0] <<" " <<data3[1] <<" ";
      } else {
        out<<enume.real()/ncalls <<" ";
        if(EstimEloc) out<<data3[0] <<" ";
      }
      out<<ovlp/ncalls <<" ";
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
                    <<LocalTimer->total("Step::Orthogonalize") <<" " 
                    <<LocalTimer->total("WalkerHandler::popControl") <<" " 
                    <<LocalTimer->total("WalkerHandler::popControl::idle") <<" " 
                    <<LocalTimer->total("WalkerHandler::loadBalance::exchange") <<" " 
                    <<max_exch_time <<" " 
                    <<LocalTimer->total("Block::TOTAL") <<" "; 
      if(prtnwalk) out<<wlkBucket->max_nbranch <<" " <<wlkBucket->max_nexch <<std::setprecision(12) <<" ";
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
    wlkBucket->nwalk_min=0;
    wlkBucket->nwalk_max=0;
    wlkBucket->max_nexch=0;
    wlkBucket->max_nbranch=0;
    ovlp=0;

    LocalTimer->reset("SubStep::Propagate");
    LocalTimer->reset("Block::EstimatorEloc");
    LocalTimer->reset("Step::PopControl"); 
    LocalTimer->reset("Step::loadBalance"); 
    LocalTimer->reset("Step::Orthogonalize"); 
    LocalTimer->reset("WalkerHandler::loadBalance::setup"); 
    LocalTimer->reset("WalkerHandler::loadBalance::exchange"); 
    LocalTimer->reset("WalkerHandler::loadBalance::resize"); 
    LocalTimer->reset("WalkerHandler::popControl");
    LocalTimer->reset("WalkerHandler::popControl::bcast");
    LocalTimer->reset("WalkerHandler::popControl::idle");
    LocalTimer->reset("Block::TOTAL"); 
    
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
    m_param.add(nwfacts,"nhist","int");
    m_param.put(cur);

    assert(nwfacts>=0);
    weight_product=ComplexType(1.0,0.0);
    for(int i=0; i<nwfacts; i++)  
      weight_factors.push(weight_product);

    app_log()<<"  BasicEstimator: Number of products in weight history: " <<nwfacts <<std::endl;

    std::transform(str1.begin(),str1.end(),str1.begin(),(int (*)(int)) tolower);
    std::transform(str2.begin(),str2.end(),str2.begin(),(int (*)(int)) tolower);
    std::transform(str3.begin(),str3.end(),str3.begin(),(int (*)(int)) tolower);
    std::transform(str4.begin(),str4.end(),str4.begin(),(int (*)(int)) tolower);
    std::transform(str5.begin(),str5.end(),str5.begin(),(int (*)(int)) tolower);
    if(str1 == "yes" || str1 == "true")
      EstimEloc = true;
    if(str2 == "no" || str2 == "false")
      timers = false;
    if(str3 == "yes" || str3 == "true")
      prtnwalk = true;
    if(str1 != "") EstimEloc_present=true; 
    if(str4 == "yes" || str4 == "true") overlap=true;
    if(str5 == "yes" || str5 == "true") diag=true;
    
    return true;
  }

  bool setup(std::vector<int>& TGdata, SPComplexSMVector *v, HamiltonianBase* ham, WavefunctionHandler* wfn,myTimer* timer, MPI_Comm heads_comm, MPI_Comm tg_comm, MPI_Comm node_comm, MPI_Comm heads_of_tg_comm) 
  {
    ncores_per_TG=TGdata[4];
    core_rank = TGdata[1]%ncores_per_TG; 
    ham0=ham;
    wfn0=wfn;
    if(!EstimEloc_present) EstimEloc = wfn0->check_initialized("Estimator");
    writer = (myComm->rank()==0);
    LocalTimer = timer;

    BlockTimer = TimerManager.createTimer("Block::EstimatorEloc");

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
    data3.resize(2);

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

    return true;
  }
  
  void setTargetWeight(RealType w0) { targetW = w0; } 

  private:

  std::vector<double> data, data2, data3;  

  int nwfacts=0;
  std::queue<ComplexType> weight_factors;
  ComplexType weight_product=ComplexType(1.0,0.0);

  ComplexType enume=0.0,edeno=0.0;
  ComplexType enume_sub=0.0,edeno_sub=0.0;
  ComplexType enume2=0.0,edeno2=0.0;
  RealType weight, weight_sub, ovlp, ovlp_sub;
  RealType targetW=1;
  int nwalk_good, nwalk, ncalls, ncalls_substep, nwalk_sub, nwalk_min, nwalk_max; 
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

  NewTimer *BlockTimer;

  int num_heads_tg;
  MPI_Comm MPI_COMM_HEAD_OF_NODES;
  MPI_Comm MPI_COMM_NODE_LOCAL;
  MPI_Comm MPI_COMM_TG_LOCAL;
  MPI_Comm MPI_COMM_TG_LOCAL_HEADS;

};
}

#endif
