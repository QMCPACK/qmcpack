//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_WAVEFUNCTIONFACTORY_H
#define QMCPLUSPLUS_AFQMC_WAVEFUNCTIONFACTORY_H

#include<iostream>
#include<vector> 
#include<map> 
#include<fstream>
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"

namespace qmcplusplus
{

namespace afqmc
{

// keep a std::map<*AFQMCInfo,SlaterDetOperations> to construct Wfns, and route all determinant operations through this object in Wfn classes


class WavefunctionFactory
{

  public:

  WavefunctionFactory(std::map<std::string,AFQMCInfo>& info):InfoMap(info)
  {
  }

  ~WavefunctionFactory()
  {
  }

  bool is_constructed(const std::string& ID) 
  {
    auto xml = xmlBlocks.find(ID);
    if(xml == xmlBlocks.end())
      APP_ABORT(" Error in WavefunctionFactory::is_constructed(string&): Missing xml block. \n");
    auto w0 = wavefunctions.find(ID);
    if( w0 == wavefunctions.end() ) 
      return false; 
    else
      return true; 
  }

  // returns a pointer to the base Wavefunction class associated with a given ID 
  Wavefunction& getWavefunction(TaskGroup_& TGprop, TaskGroup_& TGwfn, 
                                const std::string& ID, WALKER_TYPES walker_type, Hamiltonian* h, 
                                RealType cutvn=1e-6, int targetNW=1)
  {
    auto xml = xmlBlocks.find(ID);
    if(xml == xmlBlocks.end())
      APP_ABORT(" Error in WavefunctionFactory::getWavefunction(string&): Missing xml block. \n");
    auto w0 = wavefunctions.find(ID);
    if( w0 == wavefunctions.end() ) {
      auto neww = wavefunctions.insert(
                std::make_pair(ID,buildWavefunction(TGprop,TGwfn,xml->second,walker_type,
                                                    h,cutvn,targetNW)));
      if(!neww.second)
        APP_ABORT(" Error: Problems building new wavefunction in WavefunctionFactory::getWavefunction(string&). \n");
      return (neww.first)->second;  
    } else
      return w0->second;
  }

  // Use this routine to check if there is a wfn associated with a given ID
  // since getWavefunction aborts if the xml block is missing
  // returns the xmlNodePtr associated with ID
  xmlNodePtr getXML(const std::string& ID)
  {
    auto xml = xmlBlocks.find(ID);
    if(xml == xmlBlocks.end())
      return nullptr;
    else
      return xml->second;
  }

  // returns the xmlNodePtr associated with ID
  boost::multi_array<ComplexType,3>& getInitialGuess(const std::string& ID)
  {
    auto mat = initial_guess.find(ID);
    if(mat == initial_guess.end()) { 
      APP_ABORT(" Error: Missing initial guess in WavefunctionFactory. \n");
    } else 
      return mat->second;
  }

  // adds a xml block from which a Wavefunction can be built
  void push(const std::string& ID, xmlNodePtr cur)
  {
    auto xml = xmlBlocks.find(ID);
    if(xml != xmlBlocks.end())
      APP_ABORT("Error: Repeated Wavefunction block in WavefunctionFactory. Wavefunction names must be unique. \n");
    xmlBlocks.insert(std::make_pair(ID,cur));
  }

  protected:

  // reference to container of AFQMCInfo objects 
  std::map<std::string,AFQMCInfo>& InfoMap;

  // generates a new Wavefunction and returns the pointer to the base class
  Wavefunction buildWavefunction(TaskGroup_& TGprop, TaskGroup_& TGwfn, xmlNodePtr cur, 
                                 WALKER_TYPES walker_type, Hamiltonian* h, RealType cutvn, int targetNW)
  {
    std::string type;
    ParameterSet m_param;
    m_param.add(type,"filetype","std::string");
    m_param.put(cur);

    app_log()<<"\n****************************************************\n"
           <<"               Initializating Wavefunction \n"
           <<"****************************************************\n"
           <<std::endl;

    if(type == "none" || type == "ascii") {
      assert(h!=nullptr);
      return fromASCII(TGprop,TGwfn,cur,walker_type,*h,cutvn,targetNW);
    } else if(type == "hdf5")
      return fromHDF5(TGprop,TGwfn,cur,walker_type,cutvn,targetNW);
    else {
      app_error()<<"Unknown Wavefunction filetype in WavefunctionFactory::buildWavefunction(): " <<type <<std::endl;
      APP_ABORT(" Error: Unknown Wavefunction filetype in WavefunctionFactory::buildWavefunction(). \n");
    }
  }

  Wavefunction fromASCII(TaskGroup_& TGprop, TaskGroup_& TGwfn, xmlNodePtr cur, WALKER_TYPES walker_type,
                            Hamiltonian& h, RealType cutvn, int targetNW);

  Wavefunction fromHDF5(TaskGroup_& TGprop, TaskGroup_& TGwfn, xmlNodePtr cur, WALKER_TYPES walker_type, 
                            RealType cutvn, int targetNW);

  // MAM: should I store a copy rather than a pointer???
  std::map<std::string,xmlNodePtr> xmlBlocks;

  std::map<std::string,Wavefunction> wavefunctions;

  std::map<std::string,boost::multi_array<ComplexType,3>> initial_guess; 

  //std::map<AFQMCInfo,SlaterDetOperations>

};
}
}

#endif
