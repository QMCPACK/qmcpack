//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// jnkim@ornl.gov
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file qmc_common.h
 *
 * Declare static data that persists throughout simulation
 */
#ifndef QMCPLUSPLUS_GLOBAL_OBJECTS_H
#define QMCPLUSPLUS_GLOBAL_OBJECTS_H

#include <Configuration.h>

namespace qmcplusplus
{
///enumeration for main computing devices
enum {SMP=0, CUDA=1, PHI=2};

/** class to definte global variables to keep track a run
 */
struct QMCState
{
  ///true, if a run is a restart with <mcwalkerset/>
  bool is_restart;
  ///true, if density is used, e.g. MPC
  bool use_density;
  ///true, if it is a dryrun
  bool dryrun;
  ///true, if wave functions are stored for next runs
  bool save_wfs;
  ///true, if walker swap is done by async
  bool async_swap;
  ///true, print out file
  bool io_node;
  ///true, use Ewald instead of optimal breakup for the Coulomb
  bool use_ewald;
  ///int for compute_device
  int compute_device;
  ///init for <qmc/> section
  int qmc_counter;
  ///number of mpi groups
  int mpi_groups;
  /** size of memory allocated in byte per MPI
   */
  size_t memory_allocated;
  ///vacuum layer for non-periodic dimension
  double vacuum;
  ///store the name of the main eshd file name
  string master_eshd_name;

  ///constructor
  QMCState();
  ///initialize options from the command-line
  void initialize(int argc, char **argv);
  ///print command-line options
  void print_options(ostream& os);
  /** print memory increase
   * @param who the name of the class/function calling this
   * @param before memory_allocated before calling print
   */
  void print_memory_change(const string& who, size_t before);
};

///a unique QMCState during a run
extern QMCState qmc_common;
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5388 $   $Date: 2011-12-02 08:45:44 -0500 (Fri, 02 Dec 2011) $
 * $Id: Configuration.h 5388 2011-12-02 13:45:44Z jnkim $
 ***************************************************************************/
