
#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANFACTORY_H
#define QMCPLUSPLUS_AFQMC_HAMILTONIANFACTORY_H

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Hamiltonians/HamiltonianFactory_Helper.h"

namespace qmcplusplus
{
namespace afqmc
{
class HamiltonianFactory
{
public:
  HamiltonianFactory(std::map<std::string, AFQMCInfo>& info) : InfoMap(info) {}

  ~HamiltonianFactory()
  {
    // delete Hamiltonian objects
    //for (auto it = hamiltonians.begin(); it != hamiltonians.end(); ++it)
    //  delete it->second;
  }

  bool is_constructed(const std::string& ID)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      APP_ABORT(" Error in WavefunctionFactory::is_constructed(string&): Missing xml block. \n");
    auto ham = hamiltonians.find(ID);
    if (ham == hamiltonians.end())
      return false;
    else
      return true;
  }

  // returns a pointer to the base Hamiltonian class associated with a given ID
  Hamiltonian& getHamiltonian(GlobalTaskGroup& gTG, const std::string& ID)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      APP_ABORT("Error: Missing xml Block in HamiltonianFactory::getHamiltonian(string&). \n");
    auto ham = hamiltonians.find(ID);
    if (ham == hamiltonians.end())
    {
      auto newham = hamiltonians.insert(std::make_pair(ID, buildHamiltonian(gTG, xml->second)));
      if (!newham.second)
        APP_ABORT(" Error: Problems inserting new hamiltonian in HamiltonianFactory::getHamiltonian(streing&). \n");
      return (newham.first)->second;
    }
    else
      return ham->second;
  }

  // adds a xml block from which a Hamiltonian can be built
  void push(const std::string& ID, xmlNodePtr cur)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml != xmlBlocks.end())
      APP_ABORT("Error: Repeated Hamiltonian block in HamiltonianFactory. Hamiltonian names must be unique. \n");
    xmlBlocks.insert(std::make_pair(ID, cur));
  }

  xmlNodePtr getXML(const std::string& ID)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      return nullptr;
    else
      return xml->second;
  }

protected:
  // reference to container of AFQMCInfo objects
  std::map<std::string, AFQMCInfo>& InfoMap;

  // keep ownership of the TGs in the Factory
  // this way you can reuse them if necessary
  // and you don't then need to worry about the semantics of mpi3::communicator
  std::map<int, TaskGroup_> TGMap;

  // generates a new Hamiltonian and returns the pointer to the base class
  Hamiltonian buildHamiltonian(GlobalTaskGroup& gTG, xmlNodePtr cur)
  {
    std::string type("hdf5");
    ParameterSet m_param;
    m_param.add(type, "filetype");
    m_param.put(cur);

    app_log() << "\n****************************************************\n"
              << "               Initializing Hamiltonian \n"
              << "****************************************************\n"
              << std::endl;

    if (type == "hdf5")
      return fromHDF5(gTG, cur);
    else
    {
      app_error() << "Unknown Hamiltonian filetype in HamiltonianFactory::buildHamiltonian(): " << type << std::endl;
      APP_ABORT(" Error: Unknown Hamiltonian filetype in HamiltonianFactory::buildHamiltonian(). \n");
    }
    return Hamiltonian{};
  }

  Hamiltonian fromHDF5(GlobalTaskGroup& gTG, xmlNodePtr cur);

  //  Hamiltonian fromHDF5_old(GlobalTaskGroup& gTG, xmlNodePtr cur);

  TaskGroup_& getTG(GlobalTaskGroup& gTG, int nTG)
  {
    if (gTG.getTotalNodes() % nTG != 0)
      APP_ABORT("Error: number_of_TGs must divide the total number of processors. \n\n\n");
    int nnodes = gTG.getTotalNodes() / nTG;
    auto t     = TGMap.find(nnodes);
    if (t == TGMap.end())
    {
      auto p = TGMap.insert(std::make_pair(nnodes,
                                           TaskGroup_(gTG, std::string("HamiltonianTG_") + std::to_string(nnodes),
                                                      nnodes, gTG.getTotalCores())));
      if (!p.second)
        APP_ABORT(" Error: Problems creating new hamiltonian TG in HamiltonianFactory::getTG(int). \n");
      return (p.first)->second;
    }
    return t->second;
  }

  // MAM: should I store a copy rather than a pointer???
  std::map<std::string, xmlNodePtr> xmlBlocks;

  std::map<std::string, Hamiltonian> hamiltonians;

  inline HamiltonianTypes peekHamType(hdf_archive& dump)
  {
    if (dump.is_group(std::string("/Hamiltonian/KPTHC")))
      return KPTHC;
    if (dump.is_group(std::string("/Hamiltonian/THC")))
      return THC;
    if (dump.is_group(std::string("/Hamiltonian/KPFactorized")))
      return KPFactorized;
    if (dump.is_group(std::string("/Hamiltonian/DenseFactorized")))
      return RealDenseFactorized;
    if (dump.is_group(std::string("/Hamiltonian/Factorized")))
      return Factorized;
    APP_ABORT("  Error: Invalid hdf file format in peekHamType(hdf_archive). \n");
    return UNKNOWN;
  }
};
} // namespace afqmc
} // namespace qmcplusplus

#endif
