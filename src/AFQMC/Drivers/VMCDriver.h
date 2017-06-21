#ifndef QMCPLUSPLUS_AFQMC_VMCDRIVER_H
#define QMCPLUSPLUS_AFQMC_VMCDRIVER_H

#include<Message/MPIObjectBase.h>
#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Drivers/Driver.h"
#include "AFQMC/Wavefunctions/WavefunctionHandler.h"
#include "AFQMC/Propagators/PropagatorBase.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"
#include "AFQMC/Hamiltonians/HamiltonianBase.h"

namespace qmcplusplus
{

class VMCDriver: public Driver 
{

  typedef HamiltonianBase* HamPtr;
  typedef WavefunctionHandler* WfnPtr;
  typedef PropagatorBase* PropPtr;
  typedef WalkerHandlerBase* WSetPtr;
  typedef AFQMCInfo* InfoPtr;

  public:

    VMCDriver(Communicate *c):Driver(c),
        diagHam(0),diagHam_freq(10)
    {
      name = "VMC";
      project_title = "vmc";
    }

    ~VMCDriver() {}

    bool run();

    bool parse(xmlNodePtr); 

    bool setup(HamPtr,WSetPtr,PropPtr,WfnPtr);

    bool checkpoint(int,int);

    bool restart(hdf_archive&); 

    bool clear();

  protected:  

    int diagHam;
    int diagHam_freq;

};
}

#endif
