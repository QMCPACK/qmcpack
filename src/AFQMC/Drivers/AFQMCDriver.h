#ifndef QMCPLUSPLUS_AFQMC_AFQMCDRIVER_H
#define QMCPLUSPLUS_AFQMC_AFQMCDRIVER_H

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

class AFQMCDriver: public Driver 
{

  typedef HamiltonianBase* HamPtr;
  typedef WavefunctionHandler* WfnPtr;
  typedef PropagatorBase* PropPtr;
  typedef WalkerHandlerBase* WSetPtr;
  typedef AFQMCInfo* InfoPtr;

  public:

    AFQMCDriver(Communicate *c):Driver(c),
        debug(false),dShift(1.0),
        min_total_weight(0.8),accum_ovlp(false),
        diagHam(0),diagHam_freq(10),set_nWalker_target(false),
        print_timers(0),samplePeriod(-1),timer_first_call(true)
    {
      name = "AFQMC";
      project_title = "afqmc";
    }

    ~AFQMCDriver() {}

    bool run();

    bool parse(xmlNodePtr); 

    bool setup(HamPtr,WSetPtr,PropPtr,WfnPtr);

    bool checkpoint(int,int);

    bool restart(hdf_archive&); 

    bool clear();

  protected:  

    bool writeSamples();

    bool debug; 

    int print_timers;

    int samplePeriod;

    bool accum_ovlp;

    bool set_nWalker_target;

    int diagHam;
    int diagHam_freq;

//    RealType ovlp_cut;
    RealType dShift; 
    RealType min_total_weight; 
    RealType Eshift;
    RealType Etav;

    // temporary 
    LocalWalkerHandler* localwlkBucket;

    bool timer_first_call;
    std::vector<std::string> tags;
    void output_timers(std::ofstream&,int);

};
}

#endif
