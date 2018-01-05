//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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
  ///supercell in full precision
  Tensor<OHMMS_PRECISION_FULL,OHMMS_DIM> theSuperCell;

  ///store the name of the main eshd file name
  std::string master_eshd_name;

  ///constructor
  QMCState();
  ///initialize options from the command-line
  void initialize(int argc, char **argv);
  ///print command-line options
  void print_options(std::ostream& os);
  /** print memory increase
   * @param who the name of the class/function calling this
   * @param before memory_allocated before calling print
   */
  void print_memory_change(const std::string& who, size_t before);
};

///a unique QMCState during a run
extern QMCState qmc_common;
}

#endif
