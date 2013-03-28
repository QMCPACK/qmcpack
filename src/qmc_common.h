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
  struct qmc_common
  {
    ///store the name of the main eshd file name
    static string master_eshd_name;
    ///true, if density is used, e.g. MPC
    static bool use_density;
    ///true, if it is a dryrun
    static bool dryrun;
    ///true, if wave functions are stored for next runs
    static bool save_wfs;
    ///initialize options from the command-line
    static void initialize(int argc, char **argv);
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5388 $   $Date: 2011-12-02 08:45:44 -0500 (Fri, 02 Dec 2011) $
 * $Id: Configuration.h 5388 2011-12-02 13:45:44Z jnkim $ 
 ***************************************************************************/
