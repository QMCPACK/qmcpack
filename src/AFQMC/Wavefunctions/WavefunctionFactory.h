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

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"

namespace qmcplusplus
{
namespace afqmc
{
class WavefunctionFactory
{
public:
  WavefunctionFactory(std::map<std::string, AFQMCInfo>& info) : InfoMap(info) {}

  ~WavefunctionFactory() {}

  bool is_constructed(const std::string& ID)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      APP_ABORT(" Error in WavefunctionFactory::is_constructed(string&): Missing xml block. \n");
    auto w0 = wavefunctions.find(ID);
    if (w0 == wavefunctions.end())
      return false;
    else
      return true;
  }

  // returns a pointer to the base Wavefunction class associated with a given ID
  Wavefunction& getWavefunction(TaskGroup_& TGprop,
                                TaskGroup_& TGwfn,
                                const std::string& ID,
                                WALKER_TYPES walker_type,
                                Hamiltonian* h,
                                RealType cutvn = 1e-6,
                                int targetNW   = 1)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      APP_ABORT(" Error in WavefunctionFactory::getWavefunction(string&): Missing xml block. \n");
    auto w0 = wavefunctions.find(ID);
    if (w0 == wavefunctions.end())
    {
      auto neww = wavefunctions.insert(
          std::make_pair(ID, buildWavefunction(TGprop, TGwfn, xml->second, walker_type, h, cutvn, targetNW)));
      if (!neww.second)
        APP_ABORT(" Error: Problems building new wavefunction in WavefunctionFactory::getWavefunction(string&). \n");
      return (neww.first)->second;
    }
    else
      return w0->second;
  }

  // Use this routine to check if there is a wfn associated with a given ID
  // since getWavefunction aborts if the xml block is missing
  // returns the xmlNodePtr associated with ID
  xmlNodePtr getXML(const std::string& ID)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      return nullptr;
    else
      return xml->second;
  }

  // returns the xmlNodePtr associated with ID
  boost::multi::array<ComplexType, 3>& getInitialGuess(const std::string& ID)
  {
    auto mat = initial_guess.find(ID);
    if (mat == initial_guess.end())
    {
      APP_ABORT(" Error: Missing initial guess in WavefunctionFactory. \n");
    }
    return mat->second;
  }

  // adds a xml block from which a Wavefunction can be built
  void push(const std::string& ID, xmlNodePtr cur)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml != xmlBlocks.end())
      APP_ABORT("Error: Repeated Wavefunction block in WavefunctionFactory. Wavefunction names must be unique. \n");
    xmlBlocks.insert(std::make_pair(ID, cur));
  }

protected:
  // reference to container of AFQMCInfo objects
  std::map<std::string, AFQMCInfo>& InfoMap;

  // generates a new Wavefunction and returns the pointer to the base class
  Wavefunction buildWavefunction(TaskGroup_& TGprop,
                                 TaskGroup_& TGwfn,
                                 xmlNodePtr cur,
                                 WALKER_TYPES walker_type,
                                 Hamiltonian* h,
                                 RealType cutvn,
                                 int targetNW)
  {
    std::string type;
    ParameterSet m_param;
    m_param.add(type, "filetype");
    m_param.put(cur);

    app_log() << "\n****************************************************\n"
              << "               Initializing Wavefunction \n"
              << "****************************************************\n"
              << std::endl;

    if (type == "none" || type == "ascii")
    {
      assert(h != nullptr);
      return fromASCII(TGprop, TGwfn, cur, walker_type, *h, cutvn, targetNW);
    }
    else if (type == "hdf5")
      return fromHDF5(TGprop, TGwfn, cur, walker_type, *h, cutvn, targetNW);
    else
    {
      app_error() << "Unknown Wavefunction filetype in WavefunctionFactory::buildWavefunction(): " << type << std::endl;
      APP_ABORT(" Error: Unknown Wavefunction filetype in WavefunctionFactory::buildWavefunction(). \n");
    }
    return Wavefunction{};
  }

  Wavefunction fromASCII(TaskGroup_& TGprop,
                         TaskGroup_& TGwfn,
                         xmlNodePtr cur,
                         WALKER_TYPES walker_type,
                         Hamiltonian& h,
                         RealType cutvn,
                         int targetNW);

  Wavefunction fromHDF5(TaskGroup_& TGprop,
                        TaskGroup_& TGwfn,
                        xmlNodePtr cur,
                        WALKER_TYPES walker_type,
                        Hamiltonian& h,
                        RealType cutvn,
                        int targetNW);
  HamiltonianOperations getHamOps(const std::string& restart_file,
                                  WALKER_TYPES type,
                                  int NMO,
                                  int NAEA,
                                  int NAEB,
                                  std::vector<PsiT_Matrix>& PsiT,
                                  TaskGroup_& TGprop,
                                  TaskGroup_& TGwfn,
                                  RealType cutvn,
                                  RealType cutv2,
                                  int ndets_to_read,
                                  Hamiltonian& h);
  void getInitialGuess(hdf_archive& dump, std::string& name, int NMO, int NAEA, int NAEB, WALKER_TYPES walker_type);
  int getExcitation(boost::multi::array_ref<int, 1>& deti,
                    boost::multi::array_ref<int, 1>& detj,
                    std::vector<int>& excit,
                    int& perm);
  void computeVariationalEnergyPHMSD(TaskGroup_& TG,
                                     Hamiltonian& ham,
                                     boost::multi::array_ref<int, 2>& occs,
                                     std::vector<ComplexType>& coeff,
                                     int ndets,
                                     int NAEA,
                                     int NAEB,
                                     int NMO,
                                     bool recomputeCI);
  ComplexType slaterCondon0(Hamiltonian& ham, boost::multi::array_ref<int, 1>& det, int NMO);
  ComplexType slaterCondon1(Hamiltonian& ham, std::vector<int>& excit, boost::multi::array_ref<int, 1>& det, int NMO);
  ComplexType slaterCondon2(Hamiltonian& ham, std::vector<int>& excit, int NMO);


  // MAM: should I store a copy rather than a pointer???
  std::map<std::string, xmlNodePtr> xmlBlocks;

  std::map<std::string, Wavefunction> wavefunctions;

  std::map<std::string, boost::multi::array<ComplexType, 3>> initial_guess;
};
} // namespace afqmc
} // namespace qmcplusplus

#endif
