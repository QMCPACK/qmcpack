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
#include "Utilities/SimpleParser.h"
#include "Utilities/ProgressReportEngine.h"
#include "Utilities/OutputManager.h"
#include "OhmmsData/FileUtility.h"
#include "Platforms/sysutil.h"
#include "Platforms/devices.h"
#include "OhmmsApp/ProjectData.h"
#include "QMCApp/QMCMain.h"
#include "qmc_common.h"
//#include "tau/profiler.h"

void output_hardware_info(Communicate *comm, Libxml2Document &doc, xmlNodePtr root);

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
      if (c == "-disable-timers")
        TimerManager.set_timer_threshold(timer_level_none);

      // Default setting is 'timer_level_coarse'
      if (c.find("-enable-timers") < c.size())
      {
#ifndef ENABLE_TIMERS
        std::cerr << "The '-enable-timers' command line option will have no effect. This executable was built without ENABLE_TIMER set." << std::endl;
#endif
        int pos = c.find("=");
        if (pos != std::string::npos)
        {
          std::string timer_level = c.substr(pos+1);
          if (timer_level == "coarse")
          {
            TimerManager.set_timer_threshold(timer_level_coarse);
          }
          else if (timer_level == "medium")
          {
            TimerManager.set_timer_threshold(timer_level_medium);
          }
          else if (timer_level == "fine")
          {
            TimerManager.set_timer_threshold(timer_level_fine);
          }
          else
          {
            std::cerr << "Unknown timer level: " << timer_level << std::endl;
          }
        }
      }
      if (c.find("-verbosity") < c.size())
      {
        int pos = c.find("=");
        if (pos != std::string::npos)
        {
          std::string verbose_level = c.substr(pos+1);
          if (verbose_level == "low") {
            outputManager.setVerbosity(Verbosity::LOW);
          }
          else if (verbose_level == "high") {
            outputManager.setVerbosity(Verbosity::HIGH);
          }
          else if (verbose_level == "debug") {
            outputManager.setVerbosity(Verbosity::DEBUG);
          }
          else
          {
            std::cerr << "Unknown verbosity level: " << verbose_level << std::endl;
          }
        }
      }
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
    Init_CUDA();
  //safe to move on
  Communicate* qmcComm=OHMMS::Controller;
  if(inputs.size()>1)
  {
    qmcComm=new Communicate(*OHMMS::Controller,inputs.size());
    qmc_common.mpi_groups=inputs.size();
  }
  std::stringstream logname;
  int inpnum = (inputs.size() > 1) ? qmcComm->getGroupID() : 0;
  std::string myinput = inputs[qmcComm->getGroupID()];
  myinput = myinput.substr(0,myinput.size()-4);
  logname << myinput;

  if (qmcComm->rank() != 0) {
    outputManager.shutOff();
    // might need to redirect debug stream to a file per rank if debugging is enabled
  }
  if (inputs.size() > 1 && qmcComm->rank() == 0) {
    char fn[128];
    snprintf(fn, 127, "%s.g%03d.qmc",logname.str().c_str(),qmcComm->getGroupID());
    fn[127] = '\0';
    infoSummary.redirectToFile(fn);
    infoLog.redirectToSameStream(infoSummary);
    infoError.redirectToSameStream(infoSummary);
  }

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
 
  Libxml2Document timingDoc;
  timingDoc.newDoc("resources");
  output_hardware_info(qmcComm, timingDoc, timingDoc.getRoot());
  TimerManager.output_timing(qmcComm, timingDoc, timingDoc.getRoot());
  if(OHMMS::Controller->rank()==0)
  {
    timingDoc.dump(qmc->getTitle() + ".info.xml");
  }
  TimerManager.print(qmcComm);

  if(qmc)
    delete qmc;
  if(useGPU)
    Finalize_CUDA();
  OHMMS::Controller->finalize();
  return 0;
}

void output_hardware_info(Communicate *comm, Libxml2Document &doc, xmlNodePtr root)
{
  xmlNodePtr hardware = doc.addChild(root, "hardware");

  bool using_mpi = false;
#ifdef HAVE_MPI
  using_mpi = true;
  doc.addChild(hardware, "mpi_size", comm->size());
#endif
  doc.addChild(hardware, "mpi", using_mpi);

  bool using_openmp = false;
#ifdef ENABLE_OPENMP
  using_openmp = true;
  doc.addChild(hardware, "openmp_threads", omp_get_max_threads());
#endif
  doc.addChild(hardware, "openmp", using_openmp);

  bool using_gpu = false;
#ifdef QMC_CUDA
  using_gpu = true;
#endif
  doc.addChild(hardware, "gpu", using_gpu);

}
