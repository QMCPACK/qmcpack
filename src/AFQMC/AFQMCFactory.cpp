
// -*- C++ -*-
/**@file AFQMCFactory.cpp
 * @brief Top level class for AFQMC. Parses input and performs setup of classes.
 */

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<complex>
#include<tuple>
#include <queue>
#include<algorithm>
#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Configuration.h"
#include "OhmmsApp/RandomNumberControl.h"
#include <qmc_common.h>

#include "config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/AFQMCFactory.h"
#include "AFQMC/Walkers/WalkerSetFactory.hpp"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Propagators/PropagatorFactory.h"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Drivers/DriverFactory.h"

#include"AFQMC/Utilities/myTimer.h"

myTimer Timer;

namespace qmcplusplus
{

TimerList_t AFQMCTimers;
TimerNameList_t<AFQMCTimerIDs> AFQMCTimerNames =
{
  {block_timer, "Block"},
  {pseudo_energy_timer, "PseudoEnergy"},
  {energy_timer, "Energy"},
  {vHS_timer, "vHS"},
  {assemble_X_timer, "X"},
  {vbias_timer, "vbias"},
  {G_for_vbias_timer,"G_for_vbias"},
  {propagate_timer, "Propagate"},
  {back_propagate_timer, "BackPropagate"},
  {E_comm_overhead_timer, "Energy_comm_overhead"},
  {vHS_comm_overhead_timer, "vHS_comm_overhead"},
  {popcont_timer, "population_control"},
  {ortho_timer, "walker_orthogonalization"},
  {setup_timer, "setup"},
  {extra_timer, "extra"},
  {T1_t, "T1_t"},
  {T2_t, "T2_t"},
  {T3_t, "T3_t"},
  {T4_t, "T4_t"},
  {T5_t, "T5_t"},
  {T6_t, "T6_t"},
  {T7_t, "T7_t"},
  {T8_t, "T8_t"}
};

namespace afqmc
{

bool AFQMCFactory::parse(xmlNodePtr cur)
{
  if(cur == NULL)
    return false;

  xmlNodePtr curRoot=cur;
  InfoMap.clear();

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
  // Careful here, since all factories have a reference to this map
  // It must be built before any factory is used
  cur = curRoot->children;
  while (cur != NULL) {
    std::string cname((const char*)(cur->name));
    if(cname =="AFQMCInfo") {
      AFQMCInfo info;
      if(!info.parse(cur)) {
        app_error()<<"Error in AFQMCInfo::parse(xmlNodePtr)." <<std::endl;
        return false;
      }
      std::pair<std::map<std::string,AFQMCInfo>::iterator,bool> ret;
      ret = InfoMap.insert ( std::pair<std::string,AFQMCInfo>(info.name,info) );
      if (ret.second==false) {
        app_error()<<"ERROR: AFQMCInfo xml-block already defined: " <<info.name;
        return false;
      }
    }
    cur=cur->next;
  }

  // now look for non-executable blocks
  cur = curRoot->children;
  while (cur != NULL) {
    std::string cname((const char*)(cur->name));
    std::string oname("");
    OhmmsAttributeSet objAttrib;
    objAttrib.add(oname,"name");
    objAttrib.put(cur);

    if(cname =="Hamiltonian") {
      if( oname == "" ) {
        app_error()<<" Error: Missing name in xml-block:" <<cname <<std::endl;
        return false;
      }
      HamFac.push(oname,cur);
    } else if(cname == "Wavefunction") {
      if( oname == "" ) {
        app_error()<<" Error: Missing name in xml-block:" <<cname <<std::endl;
        return false;
      }
      WfnFac.push(oname,cur);
    } else if(cname == "WalkerSet") {
      if( oname == "" ) {
        app_error()<<" Error: Missing name in xml-block:" <<cname <<std::endl;
        return false;
      }
      WSetFac.push(oname,cur);
    } else if(cname == "Propagator") {
      if( oname == "" ) {
        app_error()<<" Error: Missing name in xml-block:" <<cname <<std::endl;
        return false;
      }
      PropFac.push(oname,cur);
    }
    cur = cur->next;
  }

  return true;

}

bool AFQMCFactory::execute(xmlNodePtr cur)
{
  if(cur == NULL)
    return false;

  int groupid=0; //myComm->getGroupID();
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

      // execute driver
      if(!DriverFac.executeDriver(std::string(fileroot),m_series,cur)) {
        app_error()<<"Error in DriverFactory::executeDriver::run()" <<std::endl;
        return false;
      }

      m_series++;

    }
    cur=cur->next;
  }

  return true;

}

}

}

