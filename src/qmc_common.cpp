//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <qmc_common.h>
#include <Platforms/sysutil.h>
#include "qmcpack_version.h"

namespace qmcplusplus
{
QMCState::QMCState()
{
  is_restart=false;
  use_density=false;
  dryrun=false;
  save_wfs=false;
  io_node=true;
  mpi_groups=1;
  use_ewald=false;
  qmc_counter=0;
  memory_allocated=0;
#if defined(QMC_CUDA)
  compute_device=1;
#else
  compute_device=0;
#endif
  master_eshd_name="none";
}

void QMCState::initialize(int argc, char **argv)
{
  io_node = (OHMMS::Controller->rank()==0);
  bool stopit=false;
  //going to use better option library
  int i=1;
  while(i<argc)
  {
    std::string c(argv[i]);
    if(c.find("--dryrun") < c.size())
    {
      dryrun=true;
    }
    else if(c.find("--save_wfs") < c.size())
    {
      save_wfs=(c.find("no")>=c.size());
    }
    else if(c.find("--help")< c.size())
    {
      stopit=true;
    }
    else if(c.find("--version")<c.size())
    {
      stopit=true;
    }
    else if(c.find("--vacuum")<c.size())
    {
      if(io_node)
        std::cerr << std::endl << "ERROR: command line option --vacuum has been deprecated. "
                  << "Use vacuum input tag as described in the manual." << std::endl;
      stopit=true;
    }
    else if(c.find("--noprint")<c.size())
    {//do not print Jastrow or PP
      io_node=false;
    }
    //else if(c.find("--use_ewald")<c.size())
    //{
    //  use_ewald=true;
    //}
    ++i;
  }
  if(stopit && io_node)
  {
    std::cerr << std::endl << "QMCPACK version "<< QMCPACK_VERSION_MAJOR <<"." << QMCPACK_VERSION_MINOR << "." << QMCPACK_VERSION_PATCH
        << " built on " << __DATE__ << std::endl;
#ifdef QMCPACK_GIT_BRANCH
    std::cerr << "  git branch: " << QMCPACK_GIT_BRANCH << std::endl;
    std::cerr << "  git last commit: " << QMCPACK_GIT_HASH << std::endl;
    std::cerr << "  git last commit date: " << QMCPACK_GIT_COMMIT_LAST_CHANGED << std::endl;
    std::cerr << "  git last commit subject: " << QMCPACK_GIT_COMMIT_SUBJECT << std::endl;
#endif
    std::cerr << std::endl << "Usage: qmcpack input [--dryrun --save_wfs[=no] --gpu]" << std::endl << std::endl;
  }
  if(stopit)
  {
    OHMMS::Controller->finalize();
    exit(1);
  }
}

void QMCState::print_options(std::ostream& os)
{
  os << "  Global options " << std::endl;
  if(dryrun)
    os << "  dryrun : qmc sections will be ignored." << std::endl;
  if(save_wfs)
    os << "  save_wfs=1 : save wavefunctions in hdf5. " << std::endl;
}

void QMCState::print_memory_change(const std::string& who, size_t before)
{
  before=memory_allocated-before;
  app_log() << "MEMORY increase " << (before>>20) << " MB " << who << std::endl;
}

QMCState qmc_common;

}
