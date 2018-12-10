
// -*- C++ -*-
/**@file AFQMCFactory.h
 * @brief Top level class for AFQMC. Parses input and performs setup of classes. 
 */

#ifndef QMCPLUSPLUS_AFQMCFACTORY_H
#define QMCPLUSPLUS_AFQMCFACTORY_H

//#ifdef AFQMC
#if 1>0

#include<string>
#include<vector>
#include<map>
#include <queue>
#include<algorithm>
#include<Message/MPIObjectBase.h>
#include "OhmmsApp/RandomNumberControl.h"

#include "config.h"
#include<AFQMC/Drivers/Driver.h>
#include<AFQMC/Walkers/WalkerHandlerBase.h>
#include<AFQMC/Hamiltonians/HamiltonianBase.h>
#include<AFQMC/Estimators/EstimatorHandler.h>
#include<AFQMC/Propagators/PropagatorBase.h>
#include<AFQMC/Wavefunctions/WavefunctionHandler.h>

#include "OhmmsData/libxmldefs.h"

namespace qmcplusplus
{

class AFQMCFactory: public MPIObjectBase 
{

  public:

    ///constructor
    AFQMCFactory(Communicate* c, RandomNumberControl&); 

    ///destructor
    ~AFQMCFactory() {}

    /* 
     *  Parses xml input and creates all non-executable objects. 
     *  Created objects (pointers actually) are stored in maps based on name in xml block.
     *  Executable sections (drivers) are created with objects already exiting
     *  in the maps. 
     */  
    bool parse(xmlNodePtr cur);

    /*
     *  Parses xml input and creates executable sections, using objects created during parsing. 
     */   
    bool execute(xmlNodePtr cur);

  private:

    int m_series;
    std::string project_title;

    // container of AFQMCInfo objects 
    std::map<std::string,AFQMCInfo*> InfoMap; 

    // container of walker handlers
    std::map<std::string,WalkerHandlerBase*> WalkerMap; 

    // container of hamiltonians 
    std::map<std::string,HamiltonianBase*> HamMap; 

    // container of estimators 
    std::map<std::string,EstimatorHandler*> EstimatorMap; 

    // container of propagators 
    std::map<std::string,PropagatorBase*> PropMap; 

    // container of wavefunctions
    std::map<std::string,WavefunctionHandler*> WfnMap; 

    ///random number controller
    RandomNumberControl& myRandomControl;     

};
}

#else

namespace qmcplusplus
{

class AFQMCFactory: public MPIObjectBase
{

  public:

    ///constructor
    AFQMCFactory(Communicate* c); 

    ///destructor
    ~AFQMCFactory() {} 

    /* 
     *  Parses xml input and creates all non-executable objects. 
     *  Created objects (pointers actually) are stored in maps based on name in xml block.
     *  Executable sections (drivers) are created with objects already exiting
     *  in the maps. 
     */  
    bool parse(xmlNodePtr cur);

    /*
     *  Parses xml input and creates executable sections, using objects created during parsing. 
     */   
    bool execute(xmlNodePtr cur);

};
}

#endif

#endif
