#ifndef QMCPLUSPLUS_AFQMC_AFQMCDRIVER_H
#define QMCPLUSPLUS_AFQMC_AFQMCDRIVER_H

#include "io/hdf_archive.h"
#include "mpi3/communicator.hpp"

#include "AFQMC/config.h"
#include "AFQMC/Propagators/Propagator.hpp"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Estimators/EstimatorHandler.h"
#include "Utilities/NewTimer.h"

namespace qmcplusplus
{

namespace afqmc
{

class AFQMCDriver: public AFQMCInfo
{

  public:

    AFQMCDriver(boost::mpi3::communicator& comm, AFQMCInfo& info, 
        std::string& title, int mser, int blk0, int stp0, double eshft_, xmlNodePtr cur, 
        Wavefunction& wfn_, Propagator& prpg_, EstimatorHandler& estim_):
                AFQMCInfo(info),
		globalComm(comm),
		project_title(title),
		m_series(mser),
                wfn0(wfn_),
		prop0(prpg_),
		estim0(estim_),
                block0(blk0),
		step0(stp0),
		Eshift(eshft_)
    {
      name = "AFQMCDriver";

      // read options from xml block
      parse(cur);
    }

    ~AFQMCDriver() {}

    bool run(WalkerSet&);

    bool parse(xmlNodePtr); 

    bool checkpoint(WalkerSet&,int,int);

    bool clear();

  protected:  

    boost::mpi3::communicator& globalComm;

    std::string name;

    int m_series;
    std::string project_title;

    myTimer LocalTimer;

    std::string hdf_write_restart;

    int nBlock;
    int nStep;
    int nSubstep;
    int fix_bias;

    int nCheckpoint;
    int nStabilize;
    RealType dt;
    int block0, step0;

    Wavefunction& wfn0;

    Propagator& prop0;

    EstimatorHandler& estim0;

    bool writeSamples(WalkerSet&);

    int samplePeriod;

    RealType dShift; 
    RealType Eshift;
    RealType Etav;

};

}
}

#endif
