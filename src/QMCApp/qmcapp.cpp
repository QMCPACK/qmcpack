//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Configuration.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/SimpleParser.h"
#include "OhmmsData/FileUtility.h"
#include "Platforms/sysutil.h"
#include "OhmmsApp/ProjectData.h"
#include "QMCApp/QMCMain.h"
//#include "tau/profiler.h"


#ifdef QMC_CUDA
  #include "Message/CommOperators.h"
  #include <cuda_runtime_api.h>
  #include <unistd.h>

int get_device_num()
{
  const int MAX_LEN = 200;
  int size = OHMMS::Controller->size();
  int rank = OHMMS::Controller->rank();
  
  vector<char> myname(MAX_LEN);

  gethostname(&myname[0], MAX_LEN);
  std::vector<char> host_list(MAX_LEN*size);
  for (int i=0; i<MAX_LEN; i++) 
    host_list[rank*MAX_LEN+i] = myname[i];

  OHMMS::Controller->allgather(myname, host_list, MAX_LEN);
  std::vector<std::string> hostnames;
  for (int i=0; i<size; i++) 
    hostnames.push_back(&(host_list[i*MAX_LEN]));

  string myhostname = &myname[0];
  int devnum = 0;
  for (int i=0; i<rank; i++)
    if (hostnames[i] == myhostname)
      devnum++;
  return devnum;
}

int 
get_num_appropriate_devices()
{
  
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  int num_appropriate=0;
  for (int device=0; device < deviceCount; ++device) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    if (((deviceProp.major >= 1) && (deviceProp.minor >= 3)) ||
	deviceProp.major >= 2)
      num_appropriate++;
  }
  return num_appropriate;
}

void
set_appropriate_device_num(int num)
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  int num_appropriate=0, device=0;
  for (device = 0; device < deviceCount; ++device) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    if (((deviceProp.major >= 1) && (deviceProp.minor >= 3)) ||
	deviceProp.major >= 2) {
      num_appropriate++;
      if (num_appropriate == num+1)
	cudaSetDevice (device);
    }
  }
}


void Init_CUDA(int rank, int size)
{
  int devNum = get_device_num();
  cerr << "Rank = " << rank 
       << "  My device number = " << devNum << endl;
  int num_appropriate = get_num_appropriate_devices();
  if (devNum >= num_appropriate) {
    cerr << "Not enough double-precision capable GPUs for MPI rank "
	 << rank << ".\n";
    abort();
  }
  set_appropriate_device_num (devNum);
  return;

  int numGPUs;
  cudaGetDeviceCount(&numGPUs);
  cerr << "There are " << numGPUs << " GPUs.";
  cerr << "Size = " << size << endl;
  int chunk = size/numGPUs;
  //  int device_num = rank % numGPUs;
  int device_num = rank * numGPUs / size;
  cerr << "My device number is " << device_num << endl;
  cudaSetDevice (device_num);
}
#endif


/** @file qmcapp.cpp
 *@brief a main function for QMC simulation. 
 *
 * @ingroup qmcapp
 * @brief main function for qmcapp executable.
 *
 *Actual works are done by QMCAppBase and its derived classe.
 *For other simulations, one can derive a class from QMCApps, similarly to MolecuApps.
 */
int main(int argc, char **argv) {
  ///done with the option

  //TAU_PROFILE("int main(int, char **)", " ", TAU_DEFAULT);
  //TAU_INIT(&argc, &argv);
 
  using namespace qmcplusplus;

  OHMMS::Controller->initialize(argc,argv);
  // Write out free memory on each node on Linux.

  //check the options first
  int clones=1;
  vector<string> fgroup1,fgroup2;
#ifdef QMC_CUDA
  bool useGPU = true;
#else
  bool useGPU = false;
#endif
  int i=1;
  while(i<argc)
  {
    string c(argv[i]);
    if (c.find("--gpu") < c.size()) 
      useGPU = true;
    if(c.find("clones")<c.size())
    {
      clones=atoi(argv[++i]);
    }
    else if(c.find("xml")<c.size())
    {
      fgroup1.push_back(argv[i]);
    }
    else 
    {
      ifstream fin(argv[i]);
      bool valid=true;
      do 
      {
        vector<string> words;
        getwords(words,fin);
        if(words.size())
        {
          if(words[0].find("xml")<words[0].size())
          {
            int nc=1;
            if(words.size()>1) nc=atoi(words[1].c_str());
            while(nc)
            {
              fgroup2.push_back(words[0]);--nc;
            }
          }
        }
        else
          valid=false;
      } while(valid);
    }
    ++i;
  }
  int in_files=fgroup1.size();
  vector<string> inputs(in_files*clones+fgroup2.size());
  std::copy(fgroup2.begin(),fgroup2.end(),inputs.begin());
  i=fgroup2.size();
  for(int k=0; k<in_files; ++k)
    for(int c=0; c<clones; ++c) inputs[i++]=fgroup1[k];

  if(inputs.empty())
  {
    if(OHMMS::Controller->rank()==0)
    {
      cerr << "No input file is given" << endl;
      cerr << "usage: qmcapp [--clones int] input-files " << endl;
    }
    APP_ABORT("Missing input file");
    return 1;
  }

  if (useGPU) 
#ifdef QMC_CUDA
    Init_CUDA(OHMMS::Controller->rank(),
	      OHMMS::Controller->size());
#else
  {
    cerr << "Flag \"--gpu\" was used, but QMCPACK was built without "
	 << "GPU code.\nPlease use cmake -DQMC_CUDA=1.\n";
    abort();
  }
#endif

  //safe to move on
  Communicate* qmcComm=OHMMS::Controller;
  if(inputs.size()>1)
    qmcComm=new Communicate(*OHMMS::Controller,inputs.size());

  stringstream logname;
  // logname<<getDateAndTime("%Y%m%dT%H%M");
  int inpnum = (inputs.size() > 1) ? qmcComm->getGroupID() : 0;
  string myinput = inputs[qmcComm->getGroupID()];
  myinput = myinput.substr(0,myinput.size()-4);
  //int pos = myinput.rfind(".xml");
  logname << myinput;
  OhmmsInfo Welcome(logname.str(),qmcComm->rank(),qmcComm->getGroupID(),inputs.size());

//#if defined(MPIRUN_EXTRA_ARGUMENTS)
//  //broadcast the input file name to other nodes
//  MPI_Bcast(fname.c_str(),fname.size(),MPI_CHAR,0,OHMMS::Controller->getID());
//#endif

  QMCMain *qmc=0;
  bool validInput=false;
  app_log() << "  Input file(s): ";
  for(int k=0; k<inputs.size(); ++k) app_log() << inputs[k] << " ";
  app_log() << endl;
  qmc = new QMCMain(qmcComm);
  if(inputs.size()>1)
    validInput=qmc->parse(inputs[qmcComm->getGroupID()]);
  else
    validInput=qmc->parse(inputs[0]);
  if(validInput) qmc->execute();
  if(qmc) delete qmc;

  OHMMS::Controller->finalize();
  return 0;
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
