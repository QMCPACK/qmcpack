#include <tuple>
#include <map>
#include <string>
#include <iomanip>

#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/libxmldefs.h"
#include "Configuration.h"
#include "Utilities/qmc_common.h"

#include "AFQMC/config.h"
#include "AFQMCDriver.h"
#include "AFQMC/Walkers/WalkerIO.hpp"
#include "AFQMC/Memory/buffer_managers.h"

namespace qmcplusplus
{
namespace afqmc
{
bool AFQMCDriver::run(WalkerSet& wset)
{
  std::vector<ComplexType> curData;

  RealType w0   = wset.GlobalWeight();
  int nwalk_ini = wset.GlobalPopulation();

  app_log() << "Initial weight and number of walkers: " << w0 << " " << nwalk_ini << "\n"
            << "Initial Eshift: " << Eshift << std::endl;

  // problems with using step_tot to do ortho and load balance
  double total_time = step0 * nSubstep * dt;
  int step_tot      = step0, iBlock;
  for (iBlock = block0; iBlock < nBlock; ++iBlock)
  {
    {
      ScopedTimer local_timer(AFQMCTimers[block_timer]);
      for (int iStep = 0; iStep < nStep; ++iStep, ++step_tot)
      {
        // propagate nSubstep
        prop0.Propagate(nSubstep, wset, Eshift, dt, fix_bias);
        total_time += nSubstep * dt;

        if ((step_tot + 1) % nStabilize == 0)
        {
          ScopedTimer local_timer(AFQMCTimers[ortho_timer]);
          wfn0.Orthogonalize(wset, !prop0.free_propagation());
        }

        if (total_time < weight_reset_period && !prop0.free_propagation())
          wset.resetWeights();

        {
          {
            ScopedTimer local_timer(AFQMCTimers[popcont_timer]);
            wset.popControl(curData);
          }
          estim0.accumulate_step(wset, curData);
        }

        if (total_time < 1.0)
          Eshift = estim0.getEloc_step();
        else
          Eshift += dShift * (estim0.getEloc_step() - Eshift);

        // MAM: updating here to avoid doing multiple loops with secondary allocation.
        // should do nothing after first block is finished
      }

      // checkpoint
      if (nCheckpoint > 0 && (iBlock + 1) % nCheckpoint == 0)
        if (!checkpoint(wset, iBlock, step_tot))
        {
          app_error() << " Error in AFQMCDriver::checkpoint(). \n" << std::endl;
          return false;
        }

      // write samples
      if (samplePeriod > 0 && (iBlock + 1) % samplePeriod == 0)
        if (!writeSamples(wset))
        {
          app_error() << " Error in AFQMCDriver::writeSamples(). \n" << std::endl;
          return false;
        }

      // quantities that are measured once per block
      estim0.accumulate_block(wset);
    }
    estim0.print(iBlock + 1, total_time, Eshift, wset);

    // resize stack pointers to match maximum buffer use
    update_memory_managers();
  }

  if (nCheckpoint > 0)
    checkpoint(wset, iBlock, step_tot);

  return true;
}

bool AFQMCDriver::parse(xmlNodePtr cur)
{
  if (cur == NULL)
    return false;

  // set defaults
  nBlock = 100;
  nStep = nSubstep = nStabilize = 1;
  dt                            = 0.01;
  nCheckpoint                   = -1;
  hdf_write_restart             = "";

  dShift       = 1.0;
  samplePeriod = -1;
  fix_bias     = 1;

  ParameterSet m_param;
  m_param.add(nBlock, "blocks");
  m_param.add(nStep, "steps");
  m_param.add(nSubstep, "substeps");
  m_param.add(fix_bias, "fix_bias");
  m_param.add(nStabilize, "ortho");
  m_param.add(nCheckpoint, "checkpoint");
  m_param.add(samplePeriod, "samplePeriod");
  m_param.add(weight_reset_period, "weight_reset");
  m_param.add(dt, "dt");
  m_param.add(dt, "timestep");
  m_param.add(dShift, "dshift");
  m_param.add(hdf_write_restart, "hdf_write_file");
  m_param.put(cur);

  // write all the choices here ...

  fix_bias = std::min(fix_bias, nSubstep);
  if (fix_bias > 1)
    app_log() << " Keeping the bias potential fixed for " << fix_bias << " steps. \n";

  return true;
}

// writes checkpoint file
bool AFQMCDriver::checkpoint(WalkerSet& wset, int block, int step)
{
  // hack until hdf_archive works with mpi3
  hdf_archive dump(globalComm, false);
  if (globalComm.rank() == 0)
  {
    std::string file;
    if (hdf_write_restart != std::string(""))
      file = hdf_write_restart;
    else
      file = project_title + std::string(".chk.h5");

    if (!dump.create(file))
    {
      app_error() << " Error opening checkpoint file for write. \n";
      return false;
    }

    std::vector<RealType> Rdata(2);
    Rdata[0] = Eshift;
    Rdata[1] = Eshift;

    std::vector<IndexType> Idata(2);
    Idata[0] = block;
    Idata[1] = step;

    // always write driver data and walkers
    dump.push("AFQMCDriver");
    dump.write(Idata, "DriverInts");
    dump.write(Rdata, "DriverReals");
    dump.pop();
  }

  if (!dumpToHDF5(wset, dump))
  {
    app_error() << " Problems writing checkpoint file in Driver/AFQMCDriver::checkpoint(). \n";
    return false;
  }

  if (globalComm.rank() == 0)
  {
    dump.close();
  }

  return true;
}


// writes samples
bool AFQMCDriver::writeSamples(WalkerSet& wset)
{
  // hack until hdf_archive works with mpi3
  hdf_archive dump(globalComm, false);
  if (globalComm.rank() == 0)
  {
    std::string file;
    file = project_title + std::string(".confg.h5");

    if (!dump.create(file))
    {
      app_error() << " Error opening checkpoint file for write. \n";
      return false;
    }
  }

  int nwtowrite = -1;
  if (!dumpSamplesHDF5(wset, dump, nwtowrite))
  {
    app_error() << " Problems writing checkpoint file in Driver/AFQMCDriver::writeSample(). \n";
    return false;
  }

  if (globalComm.rank() == 0)
  {
    dump.close();
  }

  return true;
}

bool AFQMCDriver::clear() { return true; }

} // namespace afqmc

} // namespace qmcplusplus
