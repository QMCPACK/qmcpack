
// -*- C++ -*-
/**@file AFQMCFactory.cpp
 * @brief Top level class for AFQMC. Parses input and performs setup of classes. 
 */

//#ifdef AFQMC
#if 1>0

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<complex>
#include<tuple>
#include <queue>
#include<algorithm>
#include<Message/MPIObjectBase.h>
#include "Message/OpenMP.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Configuration.h" 
#include "OhmmsApp/RandomNumberControl.h"
#include <qmc_common.h>

#include"AFQMC/AFQMCFactory.h"
#include "config.h"
#include<AFQMC/Walkers/WalkerHandlerBase.h>
#include<AFQMC/Walkers/DistWalkerHandler.h>
#include<AFQMC/Walkers/LocalWalkerHandler.h>
#include<AFQMC/Hamiltonians/HamiltonianBase.h>
#include<AFQMC/Estimators/EstimatorHandler.h>
#include<AFQMC/Propagators/PropagatorBase.h>
#include<AFQMC/Wavefunctions/WavefunctionHandler.h>

#include<AFQMC/Propagators/phaseless_ImpSamp_ForceBias.h>
#include<AFQMC/Propagators/VMCPropagator.h>
#include<AFQMC/Hamiltonians/SparseGeneralHamiltonian.h>
#include"AFQMC/Wavefunctions/PureSingleDeterminant.h"
#include"AFQMC/Drivers/Driver.h"
#include"AFQMC/Drivers/selectedCI.h"
#include"AFQMC/Drivers/AFQMCDriver.h"
#include"AFQMC/Drivers/BenchmarkDriver.h"
//#include"AFQMC/Drivers/AFQMCDistDriver.h"
#include"AFQMC/Drivers/VMCDriver.h"

#include"AFQMC/Utilities/myTimer.h"

// place timer here 
myTimer Timer;

namespace qmcplusplus
{


AFQMCFactory::AFQMCFactory(Communicate* c, RandomNumberControl &m):MPIObjectBase(c),myRandomControl(m),m_series(0),project_title("afqmc") 
{
}
    

bool AFQMCFactory::parse(xmlNodePtr cur)
{
/* Notes:
 *
 * 1. I don't need myComm in most classes, right now only for abort purposes. Only walker handler should need MPI
 * 2. right now hard-coding with SlaterDetWalker, later on make everything templated based on walker type.
 *
 */
  if(cur == NULL)
    return false;

  xmlNodePtr curRoot=cur;
  WalkerMap.clear();
  HamMap.clear();
  EstimatorMap.clear();
  PropMap.clear();
  WfnMap.clear();
  InfoMap.clear();

  app_log()<<" name: " <<cur->name <<std::endl;

  cur = curRoot->children;
  while (cur != NULL) {
    std::string cname((const char*)(cur->name));
    if(cname =="Project" || cname =="project") {
      OhmmsAttributeSet oAttrib;
      oAttrib.add(project_title,"id");
      oAttrib.add(project_title,"name");
      oAttrib.add(m_series,"series");
      oAttrib.put(cur);

    }
    cur=cur->next;
  }

  // first look only for AFQMCInfo
  cur = curRoot->children;
  while (cur != NULL) {
    std::string cname((const char*)(cur->name));
    if(cname =="AFQMCInfo") {
      AFQMCInfo* info = new AFQMCInfo();
      if(!info->parse(cur)) {
        app_error()<<"Error in AFQMCInfo::parse(xmlNodePtr)." <<std::endl;
        return false;
      } 
      std::pair<std::map<std::string,AFQMCInfo*>::iterator,bool> ret;
      ret = InfoMap.insert ( std::pair<std::string,AFQMCInfo*>(info->name,info) );
      if (ret.second==false) {
        app_error()<<"ERROR: AFQMCInfo xml-block already defined: " <<info->name;   
        return false;
      }      
    }
    cur=cur->next;
  }

  // now look for non-executable blocks
  cur = curRoot->children;
  while (cur != NULL) {
    std::string cname((const char*)(cur->name));

    if(cname =="Hamiltonian") {
      // building it right here since there is only 1 option
      // make a builder class that returns the pointer to the created object if necessary later
      HamiltonianBase* obj = (HamiltonianBase*) new SparseGeneralHamiltonian(myComm);
      if(!obj->parse(cur)) {
        app_error()<<"Error in SparseGeneralHamiltonian::parse(xmlNodePtr)." <<std::endl;
        return false;
      } 
      std::pair<std::map<std::string,HamiltonianBase*>::iterator,bool> ret;
      ret = HamMap.insert ( std::pair<std::string,HamiltonianBase*>(obj->name,obj) );
      if (ret.second==false) {
        app_error()<<"ERROR: HamiltonianBase xml-block already defined: " <<obj->name;
        return false;
      }
      std::string info("info0");
      OhmmsAttributeSet oAttrib;
      oAttrib.add(info,"info");
      oAttrib.put(cur);
      if(InfoMap.find(info) == InfoMap.end()) {
        app_error()<<"ERROR: Undefined info:" <<info <<"  \n";
        return false;
      }
      obj->copyInfo(*InfoMap[info]);
    } else if(cname == "Wavefunction") {
      WavefunctionHandler* obj = new WavefunctionHandler(myComm); 
      if(!obj->parse(cur)) {
        app_error()<<"Error in WavefunctionHandler::parse(xmlNodePtr)." <<std::endl;
        return false;
      }
      std::pair<std::map<std::string,WavefunctionHandler*>::iterator,bool> ret;
      ret = WfnMap.insert ( std::pair<std::string,WavefunctionHandler*>(obj->name,obj) );
      if (ret.second==false) {
        app_error()<<"ERROR: WavefunctionBase xml-block already defined: " <<obj->name;
        return false;
      }
      std::string info("info0");
      OhmmsAttributeSet oAttrib;
      oAttrib.add(info,"info");
      oAttrib.put(cur);
      if(InfoMap.find(info) == InfoMap.end()) {
        app_error()<<"ERROR: Undefined info:" <<info <<"  \n";
        return false;
      }
      obj->copyInfo(*InfoMap[info]);
    } else if(cname == "WalkerSet") {

      std::string type("distributed");
      std::string info("info0");
      OhmmsAttributeSet oAttrib;
      oAttrib.add(info,"info");
      oAttrib.add(type,"type");
      oAttrib.put(cur);

      WalkerHandlerBase* obj;
      if(type == "distributed" || type == "dist") 
        obj = (WalkerHandlerBase*) new DistWalkerHandler(myComm,RandomNumberControl::Children[0]); 
      else if(type=="local") { 
        APP_ABORT(" Error: local walker handler has been temporarily disabled. Use type=distributed.");
        obj = (WalkerHandlerBase*) new LocalWalkerHandler(myComm,RandomNumberControl::Children[0]); 
      } else {
        app_error()<<"Unknown WalkerSet type: " <<type <<std::endl;
        return false;
      }

      if(!obj->parse(cur)) {
        app_error()<<"Error in WalkerHandler::parse(xmlNodePtr)." <<std::endl;
        return false;
      }
      std::pair<std::map<std::string,WalkerHandlerBase*>::iterator,bool> ret;
      ret = WalkerMap.insert ( std::pair<std::string,WalkerHandlerBase*>(obj->name,obj) );
      if (ret.second==false) {
        app_error()<<"ERROR: WalkerHandler xml-block already defined: " <<obj->name;
        return false;
      }
      if(InfoMap.find(info) == InfoMap.end()) {
        app_error()<<"ERROR: Undefined info:" <<info <<"  \n";
        return false;
      }
      obj->copyInfo(*InfoMap[info]);
    } else if(cname == "Propagator") {

      std::string type("afqmc");
      std::string info("info0");
      OhmmsAttributeSet oAttrib;
      oAttrib.add(info,"info");
      oAttrib.add(type,"type");
      oAttrib.put(cur);

      PropagatorBase* obj;
      if(type == "afqmc") 
        obj = (PropagatorBase*) new phaseless_ImpSamp_ForceBias(myComm,RandomNumberControl::Children[0]);
//      else if(type == "vmc") 
//        obj = (PropagatorBase*) new VMCPropagator(myComm,RandomNumberControl::Children[0]);
      else {
        app_error()<<"Unknown propagator type: " <<type <<std::endl;
        return false;
      }
      if(!obj->parse(cur)) {
        app_error()<<"Error in phaseless_ImpSamp_ForceBias::parse(xmlNodePtr)." <<std::endl;
        return false;
      }
      std::pair<std::map<std::string,PropagatorBase*>::iterator,bool> ret;
      ret = PropMap.insert ( std::pair<std::string,PropagatorBase*>(obj->name,obj) );
      if (ret.second==false) {
        app_error()<<"ERROR: PropagatorBase xml-block already defined: " <<obj->name;
        return false;
      }
      if(InfoMap.find(info) == InfoMap.end()) {
        app_error()<<"ERROR: Undefined info:" <<info <<"  \n";
        return false;
      }
      obj->copyInfo(*InfoMap[info]);
    }
    cur = cur->next;
  }

    
  return true;
}

bool AFQMCFactory::execute(xmlNodePtr cur)
{
  if(cur == NULL)
    return false;

  int nproc = myComm->size();
  int nodeid = myComm->rank();
  int groupid=myComm->getGroupID();
  char fileroot[256];
  
  bool no_gtag= (qmc_common.mpi_groups==1);

  xmlNodePtr curRoot=cur;
  cur = curRoot->children;
  while (cur != NULL) {
    std::string cname((const char*)(cur->name));
    if(cname =="execute") {

      if(no_gtag) //qnproc_g == nproc)
        sprintf(fileroot,"%s.s%03d",project_title.c_str(),m_series);
      else
        sprintf(fileroot,"%s.g%03d.s%03d",project_title.c_str(),groupid,m_series);

      //set the communicator name
      myComm->setName(fileroot);

      Driver* driver; 
      // check that necessary objects exist
      std::string type("afqmc");
      std::string ham("ham0");
      std::string wfn("wfn0");
      std::string wset("wset0");
      std::string prop("prop0");
      std::string info("info0"); 
      std::string paral("no"); 
      OhmmsAttributeSet oAttrib;
      oAttrib.add(info,"info");
      oAttrib.add(prop,"prop");
      oAttrib.add(wset,"wset");
      oAttrib.add(wfn,"wfn");
      oAttrib.add(ham,"ham");
      oAttrib.add(type,"type");
      oAttrib.add(paral,"distributed");
      oAttrib.put(cur);       
      if(type == "afqmc") 
          driver = (Driver*) new AFQMCDriver(myComm);
      else if(type == "benchmark")  
        driver = (Driver*) new BenchmarkDriver(myComm);
      else if(type == "vmc")  
        driver = (Driver*) new VMCDriver(myComm);
      else if(type == "selectedCI" || type == "selectedci" || type == "selCI" || type == "selci") 
        driver = (Driver*) new selectedCI(myComm);
      else {
        app_error()<<"Unknown execute driver: " <<type <<std::endl;
        return false;
      }

      WalkerHandlerBase* wlkh0=NULL;
      HamiltonianBase* h0=NULL;
      PropagatorBase* p0=NULL;
      WavefunctionHandler* wfh0=NULL;

      if(InfoMap.find(info) == InfoMap.end()) {
        app_error()<<"ERROR: Undefined info in execute block. \n"; 
        return false;
      }
      if(WalkerMap.find(wset) != WalkerMap.end()) wlkh0 = WalkerMap[wset]; 
      if(HamMap.find(ham) != HamMap.end()) h0 = HamMap[ham];
      if(PropMap.find(prop) != PropMap.end()) p0 = PropMap[prop];
      if(WfnMap.find(wfn) != WfnMap.end()) wfh0 = WfnMap[wfn];

      driver->copyInfo(*InfoMap[info]);

      if(!driver->parse(cur)) {
        app_error()<<"Error in AFQMCDriver::parse(xmlNodePtr)." <<std::endl;
        return false;
      }

      if(!driver->setup(h0,wlkh0,p0,wfh0)) {
        app_error()<<"Error in AFQMCDriver::setup(...)." <<std::endl;
        return false;
      }      

      // execute driver
      if(!driver->run()) {
        app_error()<<"Error in AFQMCDriver::run()" <<std::endl;
        return false;
      }      

      if(!driver->clear()) {
        app_error()<<"Error in AFQMCDriver::clear()." <<std::endl;
        return false;
      }      

      m_series++;

    }
    cur=cur->next;
  }

  return true;
}

}

#else
// in case no AFQMC is compiled

#include"AFQMC/AFQMCFactory.h"
#include<iostream>

namespace qmcplusplus
{

AFQMCFactory::AFQMCFactory(Communicate* c):MPIObjectBase(c) 
{

}

bool AFQMCFactory::parse(xmlNodePtr cur)
{
  std::cerr<<"Executable not compiled with AFQMC support. \n"; 
  return false;
}

bool AFQMCFactory::execute(xmlNodePtr cur)
{
  std::cerr<<"Executable not compiled with AFQMC support. \n"; 
  return false;
}

}

#endif
