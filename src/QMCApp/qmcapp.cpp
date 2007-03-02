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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Configuration.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "OhmmsData/FileUtility.h"
#include "Utilities/SimpleParser.h"
#include "QMCApp/QMCMain.h"


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

  OHMMS::Controller->initialize(argc,argv);

  OhmmsInfo Welcome(argc,argv,OHMMS::Controller->mycontext());

  qmcplusplus::QMCMain qmc(argc,argv);

#if defined(MPIRUN_EXTRA_ARGUMENTS)
  char fname[128];
  if(OHMMS::Controller->master()) {
    sprintf(fname,"%s",argv[1]);
  }

  //broadcast the input file name to other nodes
  MPI_Bcast(fname,128,MPI_CHAR,0,OHMMS::Controller->getID());

  if(qmc.parse(fname)) {
    qmc.execute();
  }
#else
  if(argc>1) {
    string fext=getExtension(argv[1]);
    bool validInput=false;
    if(fext == "xml")
    {
      validInput=qmc.parse(argv[1]);
    }
    else
    {
      ifstream fin(argv[1]);
      vector<string> fgroup;
      bool valid=true;
      do {
        vector<string> words;
        getwords(words,fin);
        if(words.size())
          fgroup.push_back(words[0]);
        else
          valid=false;
      } while(valid);
#if defined(HAVE_MPI)
      qmc.qmcComm = new Communicate(*OHMMS::Controller,fgroup.size());
      cout << qmc.qmcComm->mycontext() << " belongs to " << qmcComm->getGroupID()
        << " using " << fgroup[qmcComm->getGroupID()] << endl;
      validInput=false;
      validInput=qmc.parse(fgroup[qmc.qmcComm->getGroupID()]);
#else
      validInput=qmc.parse(fgroup[0]);
#endif
    }
    if(validInput) qmc.execute();
  } else {
    ERRORMSG("No input file is given.")
    ERRORMSG("usage: qmcapp input-file")
  }
#endif
  LOGMSG("Bye")
  OHMMS::Controller->finalize();
  return 0;
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


