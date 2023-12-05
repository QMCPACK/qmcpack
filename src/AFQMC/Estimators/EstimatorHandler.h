#ifndef QMCPLUSPLUS_AFQMC_ESTIMATORHANDLER_H
#define QMCPLUSPLUS_AFQMC_ESTIMATORHANDLER_H

#include "Host/sysutil.h"
#include "AFQMC/config.h"

#include "AFQMC/Utilities/taskgroup.h"

#include "AFQMC/Estimators/EstimatorBase.h"
#include "AFQMC/Estimators/EnergyEstimator.h"
#include "AFQMC/Estimators/BasicEstimator.h"
#include "AFQMC/Estimators/MixedRDMEstimator.h"
#include "AFQMC/Estimators/BackPropagatedEstimator.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Hamiltonians/HamiltonianFactory.h"
#include "AFQMC/Wavefunctions/WavefunctionFactory.h"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Hamiltonians/Hamiltonian.hpp"
#include "AFQMC/Propagators/Propagator.hpp"

#include "mpi3/communicator.hpp"

namespace qmcplusplus
{
namespace afqmc
{
/* 
 * Manager class for all estimators/observables.
 * This class contains and manages a list of estimator objects.
 * An arbitrary combination of estimators can be used simultaneously
 * during a simulation, including: 
 *   1) mixed distribution estimators,  
 *   2) back propagated estimators, 
 *   3) any number of 1),2), 
 *   4) each with independent wavefunctions.
 */
class EstimatorHandler : public AFQMCInfo
{
  using EstimPtr     = std::shared_ptr<EstimatorBase>;
  using communicator = boost::mpi3::communicator;

public:
  EstimatorHandler(afqmc::TaskGroupHandler& TGgen,
                   AFQMCInfo info,
                   std::string title,
                   xmlNodePtr cur,
                   WalkerSet& wset,
                   WavefunctionFactory& WfnFac,
                   Wavefunction& wfn0,
                   Propagator& prop0,
                   WALKER_TYPES walker_type,
                   HamiltonianFactory& HamFac,
                   std::string ham0,
                   double dt,
                   bool defaultEnergyEstim = false,
                   bool impsamp            = true)
      : AFQMCInfo(info), project_title(title), hdf_output(false)
  {
    estimators.reserve(10);

    app_log() << "\n****************************************************\n"
              << "               Initializing Estimators \n"
              << "****************************************************\n"
              << std::endl;

    std::string overwrite_default_energy("no");
    xmlNodePtr curRoot  = cur;
    xmlNodePtr curBasic = NULL;
    cur                 = curRoot->children;
    while (cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if (cname == "Estimator")
      {
        std::string name("");
        OhmmsAttributeSet oAttrib;
        oAttrib.add(name, "name");
        oAttrib.put(cur);
        if (name == "basic" || name == "Basic" || name == "standard")
        {
          curBasic = cur;
        }
        else if (name == "energy")
        {
          ParameterSet m_param;
          m_param.add(overwrite_default_energy, "overwrite");
          m_param.put(cur);
        }
      }
      cur = cur->next;
    }

    estimators.emplace_back(
        static_cast<EstimPtr>(std::make_shared<BasicEstimator>(TGgen.getTG(1), info, title, curBasic, impsamp)));

    // add an EnergyEstimator if requested
    if (defaultEnergyEstim && not(overwrite_default_energy == "yes"))
      estimators.emplace_back(
          static_cast<EstimPtr>(std::make_shared<EnergyEstimator>(TGgen.getTG(1), info, nullptr, wfn0, impsamp)));


    cur = curRoot->children;
    while (cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if (cname == "Estimator")
      {
        std::string name("");
        std::string wfn_name("");
        std::string ham_name("");
        OhmmsAttributeSet oAttrib;
        oAttrib.add(name, "name");
        oAttrib.add(wfn_name, "wfn");
        oAttrib.add(ham_name, "ham");
        oAttrib.put(cur);

        if (name == "basic" || name == "Basic" || name == "standard")
        {
          // do nothing
          // first process estimators that do not need a wfn
        }
        else
        {
          // now do those that do

          Wavefunction* wfn = &wfn0;
          //Hamiltonian* ham = &ham0;
          // not sure how to do this right now
          if (wfn_name != "")
          { // wfn_name must produce a viable wfn object
            int nnodes         = 1;
            xmlNodePtr wfn_cur = WfnFac.getXML(wfn_name);
            if (wfn_cur == nullptr)
            {
              app_error() << " Error: Wavefunction named " << wfn_name << " not found. " << std::endl;
              APP_ABORT("Error: Wavefunction name not found. \n");
            }
            ParameterSet m_param;
            m_param.add(nnodes, "nnodes_per_TG");
            m_param.add(nnodes, "nnodes");
            m_param.put(wfn_cur);
            if (WfnFac.is_constructed(wfn_name))
            {
              app_warning() << " ****************************************************************************\n\n"
                            << " WARNING: Wavefunction already constructed in Estimator Handler. \n"
                            << "          Will reused available wavefunction object, ignoring hamiltonian. \n"
                            << " ****************************************************************************\n\n";
              wfn = std::addressof(
                  WfnFac.getWavefunction(TGgen.getTG(1), TGgen.getTG(nnodes), wfn_name, wfn0.getWalkerType(), nullptr));
            }
            else if (ham_name != "")
            {
              //APP_ABORT(" Estimator wfn must used default hamiltonian for execute block for now.\n");
              Hamiltonian& ham = HamFac.getHamiltonian(TGgen.gTG(), ham_name);
              wfn              = std::addressof(WfnFac.getWavefunction(TGgen.getTG(1), TGgen.getTG(nnodes), wfn_name,
                                                          wfn0.getWalkerType(), std::addressof(ham)));
            }
            else
            {
              Hamiltonian& ham = HamFac.getHamiltonian(TGgen.gTG(), ham0);
              wfn              = std::addressof(WfnFac.getWavefunction(TGgen.getTG(1), TGgen.getTG(nnodes), wfn_name,
                                                          wfn0.getWalkerType(), std::addressof(ham)));
            }
            if (wfn == nullptr)
            {
              app_error() << "Error initializing Wavefunction in DriverFactory::executeAFQMCDriver()." << std::endl;
              app_error() << "WavefunctionFactory returned nullptr, check that given Wavefunction has been defined. "
                          << std::endl;
              APP_ABORT(" Error: Problems generating wavefunction in DriverFactory::executeAFQMCDriver(). \n");
            }
          }

          if (name == "mixed_one_rdm")
          {
            estimators.emplace_back(static_cast<EstimPtr>(
                std::make_shared<MixedRDMEstimator>(TGgen.getTG(1), info, title, cur, walker_type, *wfn, impsamp)));
            hdf_output = true;
          }
          else if (name == "back_propagation")
          {
            estimators.emplace_back(static_cast<EstimPtr>(
                std::make_shared<BackPropagatedEstimator>(TGgen.getTG(1), info, title, cur, walker_type, wset, *wfn,
                                                          prop0, impsamp)));
            hdf_output = true;
          }
          else if (name == "energy")
          {
            estimators.emplace_back(
                static_cast<EstimPtr>(std::make_shared<EnergyEstimator>(TGgen.getTG(1), info, cur, *wfn, impsamp)));
          }
          else
          {
            app_log() << " Ignoring unknown estimator type: " << name << std::endl;
          }
        }
      }
      cur = cur->next;
    }


    if (TGgen.getTG(1).getGlobalRank() == 0)
    {
      //out.open(filename.c_str(),std::ios_base::app | std::ios_base::out);
      std::string filename = project_title + ".scalar.dat";
      if (hdf_output)
      {
        hdf_file = project_title + ".stat.h5";
        if (!dump.create(hdf_file))
        {
          app_log() << "Problems opening estimator hdf5 output file: " << hdf_file << std::endl;
          APP_ABORT("Problems opening estimator hdf5 output file.\n");
        }
        write_hdf_metadata(walker_type, !impsamp, dt);
      }
      out.open(filename.c_str());
      if (out.fail())
      {
        app_log() << "Problems opening estimator output file: " << filename << std::endl;
        APP_ABORT("Problems opening estimator output file. \n");
      }
      out << "# block  time  ";
      for (std::vector<EstimPtr>::iterator it = estimators.begin(); it != estimators.end(); it++)
        (*it)->tags(out);
      out << "Eshift freeMemory ";
      estimators[0]->tags_timers(out);
      out << std::endl;
    }
  }

  ~EstimatorHandler() {}

  double getEloc() { return estimators[0]->getEloc(); }

  double getEloc_step() { return estimators[0]->getEloc_step(); }

  void print(int block, double time, double Es, WalkerSet& wlks)
  {
    out << block << " " << time << " ";
    if (hdf_output)
      dump.open(hdf_file);
    for (std::vector<EstimPtr>::iterator it = estimators.begin(); it != estimators.end(); it++)
      (*it)->print(out, dump, wlks);
    out << std::setprecision(12) << Es << "  " << freemem() << " ";
    estimators[0]->print_timers(out);
    out << std::endl;
    if (hdf_output)
      dump.close();
    if ((block + 1) % 10 == 0)
      out.flush();
  }

  // 1) acumulates estimators over steps, and 2) reduces and accumulates substep estimators
  void accumulate_step(WalkerSet& wlks, std::vector<ComplexType>& curData)
  {
    for (std::vector<EstimPtr>::iterator it = estimators.begin(); it != estimators.end(); it++)
      (*it)->accumulate_step(wlks, curData);
  }

  // 1) acumulates estimators over steps, and 2) reduces and accumulates substep estimators
  void accumulate_block(WalkerSet& wlks)
  {
    for (std::vector<EstimPtr>::iterator it = estimators.begin(); it != estimators.end(); it++)
      (*it)->accumulate_block(wlks);
  }

  void write_hdf_metadata(WALKER_TYPES wlk, bool free_projection, double dt)
  {
    dump.open(hdf_file);
    dump.push("Metadata");
    dump.write(NMO, "NMO");
    dump.write(NAEA, "NAEA");
    dump.write(NAEB, "NAEB");
    int wlk_t_copy = wlk; // the actual data type of enum is implementation-defined. convert to int for file
    dump.write(wlk_t_copy, "WalkerType");
    dump.write(free_projection, "FreeProjection");
    dump.write(dt, "Timestep");
    dump.pop();
    dump.close();
  }

private:
  std::string project_title;

  std::vector<EstimPtr> estimators;
  std::vector<std::string> tags;

  std::ofstream out;
  hdf_archive dump;
  std::string hdf_file;
  bool hdf_output;
};
} // namespace afqmc
} // namespace qmcplusplus

#endif
