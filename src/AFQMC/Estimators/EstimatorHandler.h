#ifndef QMCPLUSPLUS_AFQMC_ESTIMATORHANDLER_H
#define QMCPLUSPLUS_AFQMC_ESTIMATORHANDLER_H

#include <Platforms/sysutil.h>
#include<Message/MPIObjectBase.h>
#include"AFQMC/config.h"
#include"AFQMC/Estimators/EstimatorBase.h"
#include"AFQMC/Estimators/BasicEstimator.h"
#include "AFQMC/Wavefunctions/WavefunctionHandler.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"
#include "AFQMC/Hamiltonians/HamiltonianBase.h"

namespace qmcplusplus
{

class EstimatorHandler: public MPIObjectBase, public AFQMCInfo
{
  public:

  EstimatorHandler(Communicate *c):MPIObjectBase(c),filename("energy.dat") {}

  ~EstimatorHandler() {}

  double getEloc() 
  {
    return estimators[0]->getEloc();
  } 

  double getEloc_step() 
  {
    return estimators[0]->getEloc_step();
  } 

  // meant for very cheap accumulations
  //void accumulate_substep(WalkerHandlerBase* wlks) 
  //{
    // only on basic for now
  //  estimators[0]->accumulate_substep(wlks); 
  //}

  void print(int block, double time, double Es, WalkerHandlerBase* wlks)
  {
    out<<block <<" " <<time <<" ";
    for(std::vector<EstimatorBase*>::iterator it=estimators.begin(); it!=estimators.end(); it++)
      (*it)->print(out,wlks);
    out<<std::setprecision(12) <<Es <<"  " <<freemem();
    out<<std::endl;
    if( (block+1)%10==0 ) out.flush();
  }

  // 1) acumulates estimators over steps, and 2) reduces and accumulates substep estimators 
  void accumulate_step(WalkerHandlerBase* wlks, std::vector<ComplexType>& curData) 
  {
    for(std::vector<EstimatorBase*>::iterator it=estimators.begin(); it!=estimators.end(); it++)
      (*it)->accumulate_step(wlks,curData);
  } 

  // 1) acumulates estimators over steps, and 2) reduces and accumulates substep estimators 
  void accumulate_block(WalkerHandlerBase* wlks) 
  {
    for(std::vector<EstimatorBase*>::iterator it=estimators.begin(); it!=estimators.end(); it++)
      (*it)->accumulate_block(wlks);
  } 

  void average(WalkerHandlerBase* wlks) 
  {
    for(std::vector<EstimatorBase*>::iterator it=estimators.begin(); it!=estimators.end(); it++)
      (*it)->average(wlks);
  } 

  bool parse(xmlNodePtr cur) 
  {
    if(cur==NULL) return false;

    estimators.reserve(10);
    estimators.clear();
    BasicEstimator* basic = new BasicEstimator(myComm);
    basic->copyInfo(*this);
    estimators.push_back(basic);

    xmlNodePtr curRoot=cur;
    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="Estimator") {
        OhmmsAttributeSet oAttrib;
        oAttrib.add(name,"name");
        oAttrib.put(cur);  
        if(name == "basic" || name == "Basic" || name == "standard" ) {
          if(!estimators[0]->parse(cur)) return false;
        } else {
          // create and parse new estimators 
          //->copyInfo(*this);
        }
      }
      cur = cur->next;
    } 

    return true;
  }

  // Estimator does not use TGs right now
  bool setup(std::vector<int>& TGdata, SPComplexSMVector* v, HamiltonianBase* ham0, WavefunctionHandler* wfn0, myTimer* LocalTimer, MPI_Comm heads_comm, MPI_Comm tg_comm, MPI_Comm node_comm, MPI_Comm head_tgs) 
  {

    for(std::vector<EstimatorBase*>::iterator it=estimators.begin(); it!=estimators.end(); it++)
      (*it)->setup(TGdata,v,ham0,wfn0,LocalTimer,heads_comm,tg_comm,node_comm,head_tgs);

    if(myComm->rank() == 0) {
      filename = myComm->getName()+".scalar.dat";
      //out.open(filename.c_str(),std::ios_base::app | std::ios_base::out);
      out.open(filename.c_str());
      if(out.fail()) {
        app_log()<<"Problems opening estimator output file: " <<filename <<std::endl;
        return false;
      }
      out<<"# block  time  ";
      for(std::vector<EstimatorBase*>::iterator it=estimators.begin(); it!=estimators.end(); it++)
        (*it)->tags(out);
      out<<"Eshift freeMemory ";
      out<<std::endl;
    }


    return true;
  }

  void setTargetWeight(RealType w0) { estimators[0]->setTargetWeight(w0); }
 
  private:

  std::vector<EstimatorBase*> estimators;
  std::vector<std::string> tags;

  std::string filename;

  std::ofstream out;

};
}

#endif
