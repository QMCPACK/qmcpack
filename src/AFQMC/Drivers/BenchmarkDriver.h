#ifndef QMCPLUSPLUS_AFQMC_BENCHMARKDRIVER_H
#define QMCPLUSPLUS_AFQMC_BENCHMARKDRIVER_H

#include "Message/MPIObjectBase.h"
#include "hdf/hdf_multi.h"
#include "hdf/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Wavefunctions/WavefunctionHandler.h"
#include "AFQMC/Propagators/PropagatorBase.h"
#include "AFQMC/Walkers/WalkerHandlerBase.hpp"
#include "AFQMC/Walkers/LocalWalkerHandler.h"
#include "AFQMC/Hamiltonians/HamiltonianBase.hpp"

namespace qmcplusplus
{
class BenchmarkDriver
{
  using HamPtr  = std::shared_ptr<HamiltonianBase>;
  using WSetPtr = std::shared_ptr<WalkerHandlerBase>;
  using WfnPtr  = WavefunctionHandler*;
  using PropPtr = PropagatorBase*;
  using InfoPtr = AFQMCInfo*;

public:
  BenchmarkDriver() : benchmark_list(""), maxnW(128), delnW(-1), nrepeat(5), dt(0.01)
  {
    name          = "Benchmark";
    project_title = "benchmark";
  }

  ~BenchmarkDriver() {}

  bool run();

  bool parse(xmlNodePtr);

  bool setup(HamPtr, WSetPtr, PropPtr, WfnPtr);

  bool checkpoint(int a, int b) { return true; }

  bool restart(hdf_archive& d) { return true; }

  bool clear() { return true; }

protected:
  int maxnW, delnW, nrepeat;

  RealType dt;

  std::string benchmark_list;
};
} // namespace qmcplusplus

#endif
