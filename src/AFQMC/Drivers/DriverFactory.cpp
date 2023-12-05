#include <iomanip>

#include "Configuration.h"
#include "OhmmsData/libxmldefs.h"
#include "RandomNumberControl.h"

#include "mpi3/communicator.hpp"

#include "AFQMC/Utilities/taskgroup.h"
#include "DriverFactory.h"
#include "AFQMC/Drivers/AFQMCDriver.h"
#include "AFQMC/Walkers/WalkerIO.hpp"
#include "AFQMC/Memory/buffer_managers.h"

#include "AFQMC/Walkers/WalkerSetFactory.hpp"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Propagators/PropagatorFactory.h"

#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Propagators/Propagator.hpp"
#include "AFQMC/Estimators/EstimatorHandler.h"

namespace qmcplusplus
{
namespace afqmc
{
bool DriverFactory::executeDriver(std::string title, int m_series, xmlNodePtr cur)
{
  std::string type("afqmc");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(type, "type");
  oAttrib.put(cur);

  if (type == "afqmc")
  {
    return executeAFQMCDriver(title, m_series, cur);
  }
  else if (type == "benchmark")
  {
    APP_ABORT(" FINISH \n ");
    return executeBenchmarkDriver(title, m_series, cur);
  }
  else
  {
    app_error() << "Unknown execute driver: " << type << std::endl;
    APP_ABORT(" Unknown execute driver. \n ");
    return false;
  }
}

bool DriverFactory::executeAFQMCDriver(std::string title, int m_series, xmlNodePtr cur)
{
  if (cur == NULL)
    APP_ABORT(" Error: Null xml node in DriverFactory::executeAFQMCDriver(). \n ");

  std::string ham_name("ham0");
  std::string wfn_name("wfn0");
  std::string wset_name("wset0");
  std::string prop_name("prop0");
  std::string info("info0");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(prop_name, "prop");
  oAttrib.add(wset_name, "wset");
  oAttrib.add(wfn_name, "wfn");
  oAttrib.add(ham_name, "ham");
  oAttrib.add(info, "info");
  oAttrib.put(cur);

  if (InfoMap.find(info) == InfoMap.end())
  {
    app_error() << "ERROR: Undefined info in execute block. \n";
    return false;
  }
  auto& AFinfo = InfoMap[info];
  int NMO      = AFinfo.NMO;
  int NAEB     = AFinfo.NAEB;

  std::string hdf_read_restart(""), str1("no");

  bool set_nWalker_target = false;
  double dt               = 0.01;
  int ncores_per_TG       = 1;
  int nWalkers            = 10;
  ParameterSet m_param;
  m_param.add(nWalkers, "nWalkers");
  m_param.add(ncores_per_TG, "ncores_per_TG");
  m_param.add(ncores_per_TG, "ncores");
  m_param.add(ncores_per_TG, "cores");
  m_param.add(dt, "dt");
  m_param.add(dt, "timestep");
  m_param.add(str1, "set_nWalker_to_target");
  m_param.add(str1, "set_nwalker_to_target");
  m_param.add(hdf_read_restart, "hdf_read_file");
  m_param.put(cur);

  // hard restriction for now
  bool first(false);
  if (ncores < 0)
  {
    first  = true;
    ncores = ncores_per_TG;
  }
  else if (ncores != ncores_per_TG)
    APP_ABORT(" Error: Current implementation requires the same ncores in all execution blocks. \n");

  TGHandler.setNCores(ncores);

  std::transform(str1.begin(), str1.end(), str1.begin(), (int (*)(int))tolower);
  if (str1 == "true" || str1 == "yes")
    set_nWalker_target = true;

  bool restarted = false;
  int step0      = 0;
  int block0     = 0;
  double Eshift  = 0.0;

  auto& rng = *RandomNumberControl::Children.front();

  app_log() << "\n****************************************************\n"
            << "****************************************************\n"
            << "****************************************************\n"
            << "          Beginning Driver initialization.\n"
            << "****************************************************\n"
            << "****************************************************\n"
            << "****************************************************\n"
            << std::endl;

  /*
   * Note: Hamiltonian is only needed to construct Wavefunction.
   *       If Wavefunction already exists in the factory (constructed in a previous exec block)
   *       or it is being initialized from hdf5, there is no need to build Hamiltonian.
   */

  app_log() << " Using " << ncores_per_TG << " cores per node in all TaskGroups. \n";

  hdf_archive read(gTG.Global());
  if (gTG.getGlobalRank() == 0)
  {
    if (hdf_read_restart != std::string(""))
    {
      if (read.open(hdf_read_restart, H5F_ACC_RDONLY))
      {
        // always write driver data and walkers
        try
        {
          read.push("AFQMCDriver", false);
        }
        catch (...)
        {
          return false;
        }

        std::vector<IndexType> Idata(2);
        std::vector<RealType> Rdata(2);

        if (!read.readEntry(Idata, "DriverInts"))
          return false;
        if (!read.readEntry(Rdata, "DriverReals"))
          return false;

        Eshift = Rdata[0];

        block0 = Idata[0];
        step0  = Idata[1];

        read.pop();
        restarted = true;
        read.close();
      }

      if (!restarted)
      {
        app_log() << " WARNING: Problems restarting simulation. Starting from default settings. \n";
      }
    }
  }
  gTG.Global().broadcast_value(restarted);
  if (restarted)
  {
    app_log() << " Restarted from file. Block, step: " << block0 << " " << step0 << std::endl;
    app_log() << "                      Eshift: " << Eshift << std::endl;
    gTG.Global().broadcast_value(Eshift);
    gTG.Global().broadcast_value(block0);
    gTG.Global().broadcast_value(step0);
  }

  /*
   * to do:
   *  - move all parameters related to HamiltonianOperations to wfn xml, e.g. cutvn
   *  - add logic for estimators, e.g. whether to evaluate energy, which wfn to use, etc.
   */

  // check that factories have registered objects.
  // for now, no defaults are allowed, but coming in the future
  // LATER: If missing, add a default block to the factory with the given name.
  if (WfnFac.getXML(wfn_name) == nullptr)
    APP_ABORT(" Error: Missing Wavefunction xml block. \n");
  if (PropFac.getXML(prop_name) == nullptr)
    APP_ABORT(" Error: Missing Propagator xml block. \n");
  if (WSetFac.getXML(wset_name) == nullptr)
    APP_ABORT(" Error: Missing Walker Set xml block. \n");

  // get factory parameters
  int nnodes_propg = std::max(1, get_parameter<int>(PropFac, prop_name, "nnodes", 1));
  int nnodes_wfn   = std::max(1, get_parameter<int>(WfnFac, wfn_name, "nnodes", 1));
  RealType cutvn   = get_parameter<RealType>(PropFac, prop_name, "cutoff", 1e-6);

  // setup task groups
  auto& TGprop = TGHandler.getTG(nnodes_propg);
  auto& TGwfn  = TGHandler.getTG(nnodes_wfn);

  // setting TG buffer generator here, as soon as localTG is available from any TG
  // defaults to 20MB. Read from input!!!
  std::size_t buffer_size(20);
  if (first)
    LocalTGBufferManager local_buffer(TGwfn.TG_local(), buffer_size * 1024uL * 1024uL);

  // walker set and type
  WalkerSet& wset          = WSetFac.getWalkerSet(TGHandler.getTG(1), wset_name, rng);
  WALKER_TYPES walker_type = wset.getWalkerType();

  if (not WfnFac.is_constructed(wfn_name))
  {
    // hamiltonian
    Hamiltonian& ham0 = HamFac.getHamiltonian(gTG, ham_name);

    // build wavefunction
    Wavefunction& wfn0 = WfnFac.getWavefunction(TGprop, TGwfn, wfn_name, walker_type, &ham0, cutvn, nWalkers);
  }

  // wfn builder should not use Hamiltonian pointer now
  Wavefunction& wfn0 = WfnFac.getWavefunction(TGprop, TGwfn, wfn_name, walker_type, nullptr, cutvn, nWalkers);

  // propagator
  Propagator& prop0 = PropFac.getPropagator(TGprop, prop_name, wfn0, rng);
  bool hybrid       = prop0.hybrid_propagation();

  // resize walker set
  if (restarted)
  {
    restartFromHDF5(wset, nWalkers, hdf_read_restart, read, set_nWalker_target);
    wfn0.Energy(wset);
  }
  else
  {
    auto initial_guess = WfnFac.getInitialGuess(wfn_name);
    wset.resize(
		nWalkers, 
		initial_guess[0], 
		initial_guess[1]({0, NMO}, {0, NAEB})
	);
    wfn0.Energy(wset);
    app_log() << " Energy of starting determinant: \n"
              << "  - Total energy    : " << std::setprecision(12) << wset[0].energy() << "\n"
              << "  - One-body energy : " << *wset[0].E1() << "\n"
              << "  - Coulomb energy  : " << *wset[0].EJ() << "\n"
              << "  - Exchange energy : " << *wset[0].EXX() << "\n";
    if (hybrid)
      Eshift = 0.0;
    else
      Eshift = real(wset[0].energy());
  }

  if (gTG.getGlobalRank() == 0)
    read.close();

  // is this run using importance sampling?
  bool free_proj = prop0.free_propagation();
  // if hybrid calculation, set to true
  bool addEnergyEstim = hybrid;

  // estimator setup
  EstimatorHandler estim0(TGHandler, AFinfo, title, cur, wset, WfnFac, wfn0, prop0, walker_type, HamFac, ham_name, dt,
                          addEnergyEstim, !free_proj);

  app_log() << "\n****************************************************\n"
            << "****************************************************\n"
            << "****************************************************\n"
            << "          Finished Driver initialization.\n"
            << "****************************************************\n"
            << "****************************************************\n"
            << "****************************************************\n"
            << std::endl;

  gTG.global_barrier();

  AFQMCDriver driver(gTG.Global(), AFinfo, title, m_series, block0, step0, Eshift, cur, wfn0, prop0, estim0);

  if (!driver.run(wset))
  {
    app_error() << " Problems with AFQMCDriver::run()." << std::endl;
    return false;
  }

  if (!driver.clear())
  {
    app_error() << " Problems with AFQMCDriver::clear()." << std::endl;
    return false;
  }

  return true;
}

bool DriverFactory::executeBenchmarkDriver(std::string title, int m_series, xmlNodePtr cur) { return false; }

} // namespace afqmc
} // namespace qmcplusplus
