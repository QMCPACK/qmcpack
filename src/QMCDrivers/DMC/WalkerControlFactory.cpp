//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "OhmmsData/ParameterSet.h"
#include "WalkerControlFactory.h"
#include "QMCDrivers/DMC/WalkerReconfiguration.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#if defined(HAVE_MPI)
#include "QMCDrivers/DMC/WalkerControlMPI.h"
#include "QMCDrivers/DMC/WalkerReconfigurationMPI.h"
#endif

namespace qmcplusplus
{

WalkerControlBase* createWalkerController(int nwtot, Communicate* comm, xmlNodePtr cur, bool reconfig)
{
  app_log() << "  Creating WalkerController: target  number of walkers = " << nwtot << std::endl;
  ///set of parameters
  int nmax = 0;
  std::string reconfigopt("no");
  ParameterSet m_param;
  m_param.add(nwtot, "targetWalkers");
  m_param.add(nwtot, "targetwalkers");
  m_param.add(nmax, "max_walkers");
  m_param.add(reconfigopt, "reconfiguration");
  m_param.put(cur);
  //if(nmax<0) nmax=2*nideal;
  //if(nmin<0) nmin=nideal/2;
  WalkerControlBase* wc = 0;
  int ncontexts         = comm->size();
  if (reconfigopt != "no" && reconfigopt != "runwhileincorrect")
    throw std::runtime_error("Reconfiguration is currently broken and gives incorrect results. Use dynamic "
                             "population control by setting reconfiguration=\"no\" or removing the reconfiguration "
                             "option from the DMC input section. If accessing the broken reconfiguration code path "
                             "is still desired, set reconfiguration to \"runwhileincorrect\" instead of \"yes\".");
  //bool fixw             = (reconfig || reconfigopt == "yes" || reconfigopt == "pure");
  bool fixw = (reconfig || reconfigopt == "runwhileincorrect");
  if (fixw)
  {
    int nwloc = std::max(omp_get_max_threads(), nwtot / ncontexts);
    nwtot     = nwloc * ncontexts;
  }
#if defined(HAVE_MPI)
  if (ncontexts > 1)
  {
    if (fixw)
    {
      app_log() << "  Using WalkerReconfigurationMPI for population control." << std::endl;
      wc = new WalkerReconfigurationMPI(comm);
    }
    else
    {
      app_log() << "  Using WalkerControlMPI for dynamic population control." << std::endl;
      wc = new WalkerControlMPI(comm);
    }
  }
  else
#endif
  {
    if (fixw)
    {
      app_log() << "  Using WalkerReconfiguration for population control." << std::endl;
      wc = new WalkerReconfiguration(comm);
    }
    else
    {
      app_log() << "  Using WalkerControlBase for dynamic population control." << std::endl;
      wc = new WalkerControlBase(comm);
    }
  }
  wc->set_method(fixw);
  wc->setMinMax(nwtot, nmax);
  return wc;
}


} // namespace qmcplusplus
