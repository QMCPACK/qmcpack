
// -*- C++ -*-
/**@file AFQMCFactory.cpp
 * @brief Top level class for AFQMC. Parses input and performs setup of classes.
 */

#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <complex>
#include <tuple>
#include <queue>
#include <algorithm>
#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Configuration.h"
#include "RandomNumberControl.h"
#include "Utilities/qmc_common.h"

#include "config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMCFactory.h"
#include "AFQMC/Walkers/WalkerSetFactory.hpp"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Propagators/PropagatorFactory.h"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Drivers/DriverFactory.h"

namespace qmcplusplus
{
const TimerNameList_t<AFQMCTimerIDs> AFQMCTimerNames = {{block_timer, "Block"},
                                                        {pseudo_energy_timer, "PseudoEnergy"},
                                                        {energy_timer, "Energy"},
                                                        {vHS_timer, "vHS"},
                                                        {assemble_X_timer, "X"},
                                                        {vbias_timer, "vbias"},
                                                        {G_for_vbias_timer, "G_for_vbias"},
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
                                                        {T8_t, "T8_t"}};

TimerList_t AFQMCTimers(getGlobalTimerManager(), AFQMCTimerNames, timer_level_coarse);

namespace afqmc
{
AFQMCFactory::AFQMCFactory(boost::mpi3::communicator& comm_)
    : m_series(0),
      project_title("afqmc"),
      gTG(comm_),
      TGHandler(gTG, -10),
      InfoMap(),
      HamFac(InfoMap),
      WSetFac(InfoMap),
      WfnFac(InfoMap),
      PropFac(InfoMap),
      DriverFac(gTG, TGHandler, InfoMap, WSetFac, PropFac, WfnFac, HamFac)
{
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
  // taken from src/OhmmsApp/RandomNumberControl.cpp
  int rank   = gTG.Global().rank();
  int nprocs = gTG.Global().size();
  int baseoffset;
  using uint_type = RandomNumberControl::uint_type;
  if (gTG.Global().root())
    baseoffset = static_cast<int>(static_cast<uint_type>(std::time(0)) % 1024);
  gTG.Global().broadcast_value(baseoffset);
  std::vector<uint_type> myprimes;
  RandomNumberControl::PrimeNumbers.get(baseoffset, nprocs, myprimes);
  arch::INIT(gTG.Node(), (unsigned long long int)(myprimes[rank]));
#endif
  // Global host buffers manager
  HostBufferManager host_buffer(10uL * 1024uL * 1024uL);  // setup monostate
  DeviceBufferManager dev_buffer(10uL * 1024uL * 1024uL); // setup monostate
  getGlobalTimerManager().set_timer_threshold(timer_level_coarse);
}

AFQMCFactory::~AFQMCFactory() { release_memory_managers(); }

bool AFQMCFactory::parse(xmlNodePtr cur)
{
  if (cur == NULL)
    return false;

  xmlNodePtr curRoot = cur;
  InfoMap.clear();

  cur = curRoot->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "Project" || cname == "project")
    {
      OhmmsAttributeSet oAttrib;
      oAttrib.add(project_title, "id");
      oAttrib.add(project_title, "name");
      oAttrib.add(m_series, "series");
      oAttrib.put(cur);
    }
    cur = cur->next;
  }

  // eventually read an initial buffer size from input
  // they are initialized now with 10 MBs

  // first look only for AFQMCInfo
  // Careful here, since all factories have a reference to this map
  // It must be built before any factory is used
  cur = curRoot->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "AFQMCInfo")
    {
      AFQMCInfo info;
      if (!info.parse(cur))
      {
        app_error() << "Error in AFQMCInfo::parse(xmlNodePtr)." << std::endl;
        return false;
      }
      std::pair<std::map<std::string, AFQMCInfo>::iterator, bool> ret;
      ret = InfoMap.insert(std::pair<std::string, AFQMCInfo>(info.name, info));
      if (ret.second == false)
      {
        app_error() << "ERROR: AFQMCInfo xml-block already defined: " << info.name;
        return false;
      }
    }
    cur = cur->next;
  }

  // now look for non-executable blocks
  cur = curRoot->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    std::string oname("");
    OhmmsAttributeSet objAttrib;
    objAttrib.add(oname, "name");
    objAttrib.put(cur);

    if (cname == "Hamiltonian")
    {
      if (oname == "")
      {
        app_error() << " Error: Missing name in xml-block:" << cname << std::endl;
        return false;
      }
      HamFac.push(oname, cur);
    }
    else if (cname == "Wavefunction")
    {
      if (oname == "")
      {
        app_error() << " Error: Missing name in xml-block:" << cname << std::endl;
        return false;
      }
      WfnFac.push(oname, cur);
    }
    else if (cname == "WalkerSet")
    {
      if (oname == "")
      {
        app_error() << " Error: Missing name in xml-block:" << cname << std::endl;
        return false;
      }
      WSetFac.push(oname, cur);
    }
    else if (cname == "Propagator")
    {
      if (oname == "")
      {
        app_error() << " Error: Missing name in xml-block:" << cname << std::endl;
        return false;
      }
      PropFac.push(oname, cur);
    }
    cur = cur->next;
  }

  return true;
}

bool AFQMCFactory::execute(xmlNodePtr cur)
{
  if (cur == nullptr)
    return false;

  int groupid = 0; //myComm->getGroupID();
  std::array<char, 256> fileroot;

  bool no_gtag = (qmc_common.mpi_groups == 1);

  xmlNodePtr curRoot = cur;
  cur                = curRoot->children;
  while (cur != nullptr)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "execute")
    {
      int length{0};
      if (no_gtag) //qnproc_g == nproc)
        length = std::snprintf(fileroot.data(), fileroot.size(), "%s.s%03d", project_title.c_str(), m_series);
      else
        length =
            std::snprintf(fileroot.data(), fileroot.size(), "%s.g%03d.s%03d", project_title.c_str(), groupid, m_series);

      if (length < 0)
        throw std::runtime_error("Error generating fileroot");

      // execute driver
      if (!DriverFac.executeDriver(std::string(fileroot.data(), length), m_series, cur))
      {
        app_error() << "Error in DriverFactory::executeDriver::run()" << std::endl;
        return false;
      }

      m_series++;
    }
    cur = cur->next;
  }

  return true;
}

} // namespace afqmc

} // namespace qmcplusplus
