//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <qmc_common.h>
#include <Platforms/sysutil.h>

namespace qmcplusplus
{
QMCState::QMCState()
{
  is_restart=false;
  use_density=false;
  dryrun=false;
  save_wfs=false;
  async_swap=false;
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
  vacuum=1.0;
  master_eshd_name="none";
}

void QMCState::initialize(int argc, char **argv)
{
  io_node= (OHMMS::Controller->rank()==0);
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
    else if(c.find("--async_swap") < c.size())
    {
      async_swap=(c.find("no")>=c.size());
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
      vacuum=atof(argv[++i]);
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
  if(stopit)
  {
    std::cerr << std::endl << "QMCPACK version "<< QMCPLUSPLUS_VERSION_MAJOR <<"." << QMCPLUSPLUS_VERSION_MINOR << "." << QMCPLUSPLUS_VERSION_PATCH
        << " subversion " << QMCPLUSPLUS_BRANCH
        << " build on " << getDateAndTime("%Y%m%d_%H%M") << std::endl;
    std::cerr << "Usage: qmcpack input [--dryrun --save_wfs[=no] --async_swap[=no] --gpu]" << std::endl << std::endl;
  }
}

void QMCState::print_options(std::ostream& os)
{
  os << "  Global options " << std::endl;
  if(dryrun)
    os << "  dryrun : qmc sections will be ignored." << std::endl;
  if(save_wfs)
    os << "  save_wfs=1 : save wavefunctions in hdf5. " << std::endl;
  if(async_swap)
    os << "  async_swap=1 : using async isend/irecv for walker swaps " << std::endl;
  else
    os << "  async_swap=0 : using blocking send/recv for walker swaps " << std::endl;
}

void QMCState::print_memory_change(const std::string& who, size_t before)
{
  before=memory_allocated-before;
  app_log() << "MEMORY increase " << (before>>20) << " MB " << who << std::endl;
}

QMCState qmc_common;

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5388 $   $Date: 2011-12-02 08:45:44 -0500 (Fri, 02 Dec 2011) $
 * $Id: Configuration.h 5388 2011-12-02 13:45:44Z jnkim $
 ***************************************************************************/
