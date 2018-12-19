#ifndef QMCPLUSPLUS_AFQMC_ESTIMATORHANDLER_H
#define QMCPLUSPLUS_AFQMC_ESTIMATORHANDLER_H

#include <Platforms/sysutil.h>
#include"AFQMC/config.h"

#include "AFQMC/Utilities/taskgroup.h"

#include "AFQMC/Estimators/EstimatorBase.h"
#include "AFQMC/Estimators/EnergyEstimator.h"
#include "AFQMC/Estimators/BasicEstimator.h"
//#include "AFQMC/Estimators/OneRdmEstimator.h"
//#include "AFQMC/Estimators/BackPropagatedEstimator.h"
//#include "AFQMC/Estimators/WalkerDMEstimator.h"

#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"

#include "mpi3/communicator.hpp"

namespace qmcplusplus
{

namespace afqmc
{

class EstimatorHandler: public AFQMCInfo
{

  using EstimPtr = std::shared_ptr<EstimatorBase>; 
  using communicator = boost::mpi3::communicator;

  public:

  EstimatorHandler(afqmc::TaskGroupHandler& TGgen, AFQMCInfo info, std::string title, xmlNodePtr cur, 
        WavefunctionFactory& WfnFac, 
        Wavefunction& wfn0, 
        Hamiltonian* ham0, 
        bool defaultEnergyEstim=false, 
        bool impsamp=true):
            AFQMCInfo(info),
            project_title(title) 
  {
    estimators.reserve(10);

    xmlNodePtr curRoot = cur; 
    xmlNodePtr curBasic = NULL; 
    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="Estimator") {
        std::string name("");
        OhmmsAttributeSet oAttrib;
        oAttrib.add(name,"name");
        oAttrib.put(cur);
        if(name == "basic" || name == "Basic" || name == "standard" ) {
          curBasic = cur;
          break;   
        } 
      }
      cur = cur->next;
    }

    estimators.emplace_back(static_cast<EstimPtr>(std::make_shared<BasicEstimator>(TGgen.getTG(1),info,title,curBasic,impsamp)));

    // add an EnergyEstimator if requested
    if(defaultEnergyEstim) 
      estimators.emplace_back(static_cast<EstimPtr>(std::make_shared<EnergyEstimator>(TGgen.getTG(1),info,nullptr,wfn0,impsamp)));


    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="Estimator") {
        std::string name("");
        std::string wfn_name("");
        OhmmsAttributeSet oAttrib;
        oAttrib.add(name,"name");
        oAttrib.add(wfn_name,"wfn");
        oAttrib.put(cur);

        // first process estimators that do not need a wfn
        if (name == "walker_density_matrix") {
//          estimators.emplace_back(static_cast<EstimPtr>(std::make_shared<WalkerDMEstimator>(TGgen.getTG(1),info,title,cur)));
        } else {
        // now do those that do

          Wavefunction* wfn = &wfn0;
/*
// not sure how to do this right now
          if(wfn_name != "") { // wfn_name must produce a viable wfn object
            int nnodes=1;
            xmlNodePtr wfn_cur = WfnFac.getXML(wfn_name);
            if(wfn_cur == nullptr) {
              app_error()<<" Error: Wavefunction named " <<wfn_name <<" not found. " <<std::endl;
              APP_ABORT("Error: Wavefunction name not found. \n");
            }
            ParameterSet m_param;
            m_param.add(nnodes,"nnodes_per_TG","int");
            m_param.add(nnodes,"nnodes","int");
            m_param.put(wfn_cur);
            wfn = WfnFac.getWavefunction(TGgen.getTG(nnodes),wfn_name,ham0);
            if(wfn == nullptr) {
              app_error()<<"Error initializing Wavefunction in DriverFactory::executeAFQMCDriver()." <<std::endl;
              app_error()<<"WavefunctionFactory returned nullptr, check that given Wavefunction has been defined. " <<std::endl;
              APP_ABORT(" Error: Problems generating wavefunction in DriverFactory::executeAFQMCDriver(). \n");
            }
          }
*/

          if (name == "reduced_density_matrix") {
//            estimators.emplace_back(static_cast<EstimPtr>(std::make_shared<OneRdmEstimator>(TGgen.getTG(1),info,title,cur,*wfn)));
          } else if (name == "back_propagation") {
//            estimators.emplace_back(static_cast<EstimPtr>(std::make_shared<BackPropagatedEstimator>(TGgen.getTG(1),info,title,cur,*wfn)));
          } else if (name == "energy") {
            estimators.emplace_back(static_cast<EstimPtr>(std::make_shared<EnergyEstimator>(TGgen.getTG(1),info,cur,*wfn,impsamp)));
          } else {
            app_log()<<" Ignoring unknown estimator type: " <<name <<std::endl;
          }
        }
      }
      cur = cur->next;
    }


    if(TGgen.getTG(1).getGlobalRank() == 0) {
      //out.open(filename.c_str(),std::ios_base::app | std::ios_base::out);
      std::string filename = project_title+".scalar.dat";
      out.open(filename.c_str());
      if(out.fail()) {
        app_log()<<"Problems opening estimator output file: " <<filename <<std::endl;
        APP_ABORT("Problems opening estimator output file. \n");
      }
      out<<"# block  time  ";
      for(std::vector<EstimPtr>::iterator it=estimators.begin(); it!=estimators.end(); it++)
        (*it)->tags(out);
      out<<"Eshift freeMemory ";
      estimators[0]->tags_timers(out);
      out<<std::endl;
    }
  

  }

  ~EstimatorHandler() {}

  double getEloc() 
  {
    return estimators[0]->getEloc();
  } 

  double getEloc_step() 
  {
    return estimators[0]->getEloc_step();
  } 

  void print(int block, double time, double Es, WalkerSet& wlks)
  {
    out<<block <<" " <<time <<" ";
    for(std::vector<EstimPtr>::iterator it=estimators.begin(); it!=estimators.end(); it++)
      (*it)->print(out,wlks);
    out<<std::setprecision(12) <<Es <<"  " <<freemem() <<" ";
    estimators[0]->print_timers(out);
    out<<std::endl;
    if( (block+1)%10==0 ) out.flush();
  }

  // 1) acumulates estimators over steps, and 2) reduces and accumulates substep estimators 
  void accumulate_step(WalkerSet& wlks, std::vector<ComplexType>& curData) 
  {
    for(std::vector<EstimPtr>::iterator it=estimators.begin(); it!=estimators.end(); it++)
      (*it)->accumulate_step(wlks,curData);
  } 

  // 1) acumulates estimators over steps, and 2) reduces and accumulates substep estimators 
  void accumulate_block(WalkerSet& wlks) 
  {
    for(std::vector<EstimPtr>::iterator it=estimators.begin(); it!=estimators.end(); it++)
      (*it)->accumulate_block(wlks);
  } 

  private:

  std::string project_title;

  std::vector<EstimPtr> estimators;
  std::vector<std::string> tags;

  std::ofstream out;

};
}
}

#endif
