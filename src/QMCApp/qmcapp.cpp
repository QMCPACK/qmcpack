//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Configuration.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/SimpleParser.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/FileUtility.h"
#include "Platforms/sysutil.h"
#include "Platforms/devices.h"
#include "OhmmsApp/ProjectData.h"
#include "QMCApp/QMCMain.h"
#include "qmc_common.h"
//#include "tau/profiler.h"


/** @file qmcapp.cpp
 *@brief a main function for QMC simulation.
 *
 * @ingroup qmcapp
 * @brief main function for qmcapp executable.
 *
 *Actual works are done by QMCAppBase and its derived classe.
 *For other simulations, one can derive a class from QMCApps, similarly to MolecuApps.
 */
int main(int argc, char **argv)
{
  //TAU_PROFILE("int main(int, char **)", " ", TAU_DEFAULT);
  //TAU_INIT(&argc, &argv);
  using namespace qmcplusplus;
  //qmc_common  and MPI is initialized
  OHMMS::Controller->initialize(argc,argv);
  int clones=1;
  bool useGPU=(qmc_common.compute_device == 1);
  std::vector<std::string> fgroup1,fgroup2;
  int i=1;
  while(i<argc)
  {
    std::string c(argv[i]);
    if(c[0]=='-')
    {
      if (c.find("gpu") < c.size())
        useGPU = true;
      if(c.find("clones")<c.size())
        clones=atoi(argv[++i]);
      if (c == "-debug")
        ReportEngine::enableOutput();
    }
    else
    {
      if(c.find("xml")<c.size())
        fgroup1.push_back(argv[i]);
      else
      {
        std::ifstream fin(argv[i],std::ifstream::in);
        bool valid=!fin.fail();
        while(valid)
        {
          std::vector<std::string> words;
          getwords(words,fin);
          if(words.size())
          {
            if(words[0].find("xml")<words[0].size())
            {
              int nc=1;
              if(words.size()>1)
                nc=atoi(words[1].c_str());
              while(nc)
              {
                fgroup2.push_back(words[0]);
                --nc;
              }
            }
          }
          else
            valid=false;
        }
      }
    }
    ++i;
  }
  int in_files=fgroup1.size();
  std::vector<std::string> inputs(in_files*clones+fgroup2.size());
  copy(fgroup2.begin(),fgroup2.end(),inputs.begin());
  i=fgroup2.size();
  for(int k=0; k<in_files; ++k)
    for(int c=0; c<clones; ++c)
      inputs[i++]=fgroup1[k];
  if(inputs.empty())
  {
    if(OHMMS::Controller->rank()==0)
    {
      std::cerr << "No input file is given." << std::endl;
      std::cerr << "Usage: qmcpack <input-files> " << std::endl;
    }
    OHMMS::Controller->finalize();
    return 1;
  }
  if (useGPU)
    Init_CUDA(OHMMS::Controller->rank(), OHMMS::Controller->size());
  //safe to move on
  Communicate* qmcComm=OHMMS::Controller;
  if(inputs.size()>1)
    qmcComm=new Communicate(*OHMMS::Controller,inputs.size());
  std::stringstream logname;
  int inpnum = (inputs.size() > 1) ? qmcComm->getGroupID() : 0;
  std::string myinput = inputs[qmcComm->getGroupID()];
  myinput = myinput.substr(0,myinput.size()-4);
  logname << myinput;
  OhmmsInfo Welcome(logname.str(),qmcComm->rank(),qmcComm->getGroupID(),inputs.size());
//#if defined(MPIRUN_EXTRA_ARGUMENTS)
//  //broadcast the input file name to other nodes
//  MPI_Bcast(fname.c_str(),fname.size(),MPI_CHAR,0,OHMMS::Controller->getID());
//#endif
  QMCMain *qmc=0;
  bool validInput=false;
  app_log() << "  Input file(s): ";
  for(int k=0; k<inputs.size(); ++k)
    app_log() << inputs[k] << " ";
  app_log() << std::endl;
  qmc = new QMCMain(qmcComm);
  if(inputs.size()>1)
    validInput=qmc->parse(inputs[qmcComm->getGroupID()]);
  else
    validInput=qmc->parse(inputs[0]);
  if(validInput)
    qmc->execute();
 
  TimerManager.output_timing(qmcComm, qmc->getTitle());
  TimerManager.print(qmcComm);

  if(qmc)
    delete qmc;
  if(useGPU)
    Finalize_CUDA();
  OHMMS::Controller->finalize();
  return 0;
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
