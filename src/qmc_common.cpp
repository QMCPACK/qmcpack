//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// jnkim@ornl.gov
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
    string c(argv[i]);
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
    cerr<<endl << "QMCPACK version "<< QMCPLUSPLUS_VERSION_MAJOR <<"." << QMCPLUSPLUS_VERSION_MINOR << "." << QMCPLUSPLUS_VERSION_PATCH
        << " subversion " << QMCPLUSPLUS_BRANCH
        << " build on " << getDateAndTime("%Y%m%d_%H%M") << endl;
    cerr << "Usage: qmcapp input [--dryrun --save_wfs[=no] --async_swap[=no] --gpu]" << endl << endl;
  }
}

void QMCState::print_options(ostream& os)
{
  os << "  Global options " << endl;
  if(dryrun)
    os << "  dryrun : qmc sections will be ignored." << endl;
  if(save_wfs)
    os << "  save_wfs=1 : save wavefunctions in hdf5. " << endl;
  if(async_swap)
    os << "  async_swap=1 : using async isend/irecv for walker swaps " << endl;
  else
    os << "  async_swap=0 : using blocking send/recv for walker swaps " << endl;
}

void QMCState::print_memory_change(const string& who, size_t before)
{
  before=memory_allocated-before;
  app_log() << "MEMORY increase " << (before>>20) << " MB " << who << endl;
}

QMCState qmc_common;

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5388 $   $Date: 2011-12-02 08:45:44 -0500 (Fri, 02 Dec 2011) $
 * $Id: Configuration.h 5388 2011-12-02 13:45:44Z jnkim $
 ***************************************************************************/
