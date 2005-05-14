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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "QMCApp/MO2Grid3D.h"

/**file main.cpp
 *@brief a main function for QMC simulation. 
 *Actual works are done by QMCApps and its derived classe.
 *For other simulations, one can derive a class from QMCApps, similarly to MolecuApps.
 *
 *Only requirements are the constructor and init function to initialize.
 */
int main(int argc, char **argv) {

  OHMMS::Controller->initialize(argc,argv);

  OhmmsInfo welcome(argc,argv,OHMMS::Controller->mycontext());

  ohmmsqmc::MO2Grid3D qmc(argc,argv);

  if(argc>1) {
    if(qmc.parse(argv[1])) {
      qmc.execute();
      qmc.saveXml();
    }
    //xmlFreeDoc(m_doc);
  } else {
    ERRORMSG("No input file is given.")
    ERRORMSG("usage: mo2grid3d input-file")
  }
  LOGMSG("Bye")

  OHMMS::Controller->finalize();
  return 0;
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


