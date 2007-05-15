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

  if(argc<=1)
  {
    ERRORMSG("No input file is given.")
    ERRORMSG("usage: qmcapp input-file")
    return 1;
  }

  qmcplusplus::QMCMain qmc(argc,argv);

  string fname=argv[1];

#if defined(MPIRUN_EXTRA_ARGUMENTS)
  //broadcast the input file name to other nodes
  MPI_Bcast(fname.c_str(),fname.size(),MPI_CHAR,0,OHMMS::Controller->getID());
#endif

  string fext=getExtension(fname);
  bool validInput=false;

  if(fext == "xml")
  {
    qmc.qmcComm = OHMMS::Controller;
    validInput=qmc.parse(fname);
  }
  else
  {
    ifstream fin(fname.c_str());
    vector<string> fgroup;
    bool valid=true;
    do 
    {
      vector<string> words;
      getwords(words,fin);
      if(words.size())
        fgroup.push_back(words[0]);
      else
        valid=false;
    } while(valid);

    if(OHMMS::Controller->ncontexts()==1)
    {
      qmc.qmcComm = OHMMS::Controller;
      validInput=qmc.parse(fgroup[0]);
    }
    else
    {
      qmc.qmcComm = new Communicate(*OHMMS::Controller,fgroup.size());
      validInput=qmc.parse(fgroup[qmc.qmcComm->getGroupID()]);
    }
  }

  if(validInput) qmc.execute();

  LOGMSG("Bye")
  OHMMS::Controller->finalize();
  return 0;
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
