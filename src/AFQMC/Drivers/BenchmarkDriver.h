#ifndef QMCPLUSPLUS_AFQMC_BENCHMARKDRIVER_H
#define QMCPLUSPLUS_AFQMC_BENCHMARKDRIVER_H

#include<Message/MPIObjectBase.h>
#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Drivers/Driver.h"
#include "AFQMC/Wavefunctions/WavefunctionHandler.h"
#include "AFQMC/Propagators/PropagatorBase.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"
#include "AFQMC/Walkers/LocalWalkerHandler.h"
#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "Utilities/NewTimer.h"

namespace qmcplusplus
{

class BenchmarkDriver: public Driver 
{

  typedef HamiltonianBase* HamPtr;
  typedef WavefunctionHandler* WfnPtr;
  typedef PropagatorBase* PropPtr;
  typedef WalkerHandlerBase* WSetPtr;
  typedef AFQMCInfo* InfoPtr;

  public:

    BenchmarkDriver(Communicate *c):Driver(c),benchmark_list("")
            ,maxnW(128),delnW(-1),nrepeat(5),dt(0.01)
    {
      name = "Benchmark";
      project_title = "benchmark";
    }

    ~BenchmarkDriver() {}

    bool run();

    bool parse(xmlNodePtr); 

    bool setup(HamPtr,WSetPtr,PropPtr,WfnPtr);

    bool checkpoint(int a,int b) { return true; }

    bool restart(hdf_archive& d) { return true; }

    bool clear() { return true; }

  protected:  

    int maxnW,delnW,nrepeat; 

    RealType dt;

    std::string benchmark_list;

};
}

#endif
