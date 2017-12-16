#ifndef QMCPLUSPLUS_AFQMC_DRIVER_H
#define QMCPLUSPLUS_AFQMC_DRIVER_H

#include<Message/MPIObjectBase.h>
#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Wavefunctions/WavefunctionHandler.h"
#include "AFQMC/Propagators/PropagatorBase.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"
#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Estimators/EstimatorHandler.h"
#include "AFQMC/Utilities/taskgroup.h"

namespace qmcplusplus
{

class Driver: public MPIObjectBase, public AFQMCInfo
{

  public:

    typedef HamiltonianBase* HamPtr;
    typedef WavefunctionHandler* WfnPtr;
    typedef PropagatorBase* PropPtr;
    typedef WalkerHandlerBase* WSetPtr;
    typedef AFQMCInfo* InfoPtr;

    Driver(Communicate *c):MPIObjectBase(c),TG(c,"DriverTG"),
        nBlock(100),nStep(1),nSubstep(1),
        nStabalize(1),nPopulationControl(1),nloadBalance(1),dt(0.01),name(""),
        nCheckpoint(-1),block0(0),step0(0),restarted(false),
        hdf_write_restart(""),hdf_read_restart(""),nWalkers(5),
        hdf_write_tag(""),hdf_read_tag(""),project_title(""),m_series(0),
        ncores_per_TG(1)
    {}

    ~Driver() {}

    virtual bool run()=0;

    void setTitle(std::string& title, int cnt) {
      project_title=title;
      m_series=cnt;
    }

    virtual bool parse(xmlNodePtr)=0; 

    virtual bool setup(HamPtr,WSetPtr,PropPtr,WfnPtr)=0;

    virtual bool checkpoint(int,int)=0;

    virtual bool restart(hdf_archive&)=0; 

    virtual bool clear()=0;

    std::string name;

  protected:  

    int m_series;
    std::string project_title;  

    myTimer LocalTimer;

    std::string hdf_read_restart;
    std::string hdf_write_restart;
    std::string hdf_read_tag;
    std::string hdf_write_tag;

    bool restarted;

    int nBlock;
    int nStep;
    int nSubstep;

    afqmc::TaskGroup TG;

    int ncores_per_TG;

    int nWalkers;

    int nCheckpoint;
    int nStabalize;
    int nPopulationControl;
    int nloadBalance;
    RealType dt;
    int block0, step0; 

    HamPtr ham0; 

    WfnPtr wfn0;

    WSetPtr wlkBucket;

    PropPtr prop0;

    EstimatorHandler* estim0;
   
    SPComplexSMVector CommBuffer; 

    MPI_Comm MPI_COMM_HEAD_OF_NODES;
    MPI_Comm MPI_COMM_NODE_LOCAL;
    MPI_Comm MPI_COMM_TG_LOCAL;
    MPI_Comm MPI_COMM_TG_LOCAL_HEADS;

};
}

#endif
