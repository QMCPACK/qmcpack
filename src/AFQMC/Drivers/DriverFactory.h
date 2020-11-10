
#ifndef QMCPLUSPLUS_DRIVERFACTORY_H
#define QMCPLUSPLUS_DRIVERFACTORY_H

#include "OhmmsData/libxmldefs.h"
#include "Message/MPIObjectBase.h"

#include "mpi3/communicator.hpp"

#include "AFQMC/Utilities/taskgroup.h"

#include "AFQMC/Walkers/WalkerSetFactory.hpp"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Propagators/PropagatorFactory.h"

namespace qmcplusplus
{
namespace afqmc
{
class DriverFactory
{
  using communicator = boost::mpi3::communicator;

public:
  DriverFactory(GlobalTaskGroup& gtg_,
                TaskGroupHandler& tghandler_,
                std::map<std::string, AFQMCInfo>& info,
                WalkerSetFactory& wsetfac_,
                PropagatorFactory& pfac_,
                WavefunctionFactory& wfnfac_,
                HamiltonianFactory& hfac)
      : ncores(-10),
        gTG(gtg_),
        TGHandler(tghandler_),
        InfoMap(info),
        WSetFac(wsetfac_),
        PropFac(pfac_),
        HamFac(hfac),
        WfnFac(wfnfac_)
  {}

  ~DriverFactory() {}


  bool executeDriver(std::string title, int m_series, xmlNodePtr cur);

private:
  bool executeAFQMCDriver(std::string title, int m_series, xmlNodePtr cur);
  bool executeBenchmarkDriver(std::string title, int m_series, xmlNodePtr cur);

  int ncores;

  // global TG from which all TGs are constructed
  GlobalTaskGroup& gTG;

  TaskGroupHandler& TGHandler;

  // container of AFQMCInfo objects
  std::map<std::string, AFQMCInfo>& InfoMap;

  // WalkerHandler factory
  WalkerSetFactory& WSetFac;

  // Propagator factory
  PropagatorFactory& PropFac;

  // Hamiltonian factory
  HamiltonianFactory& HamFac;

  // Wavefunction factory
  WavefunctionFactory& WfnFac;
};

} // namespace afqmc

} // namespace qmcplusplus

#endif
