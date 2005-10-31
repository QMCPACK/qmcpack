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
#include "QMC/MolecuApps.h"

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

  qmcplusplus::MolecuApps qmc(argc,argv);

  if(argc>1) {
    // build an XML tree from a the file;
    xmlDocPtr m_doc = xmlParseFile(argv[1]);
    if (m_doc == NULL) {
      ERRORMSG("File " << argv[1] << " is invalid")
      xmlFreeDoc(m_doc);
      return 1;
    }    
    // Check the document is of the right kind
    xmlNodePtr cur = xmlDocGetRootElement(m_doc);
    if (cur == NULL) {
      ERRORMSG("Empty document");
      xmlFreeDoc(m_doc);
      return 1;
    }
    qmc.run(cur);
  } else {
    WARNMSG("No argument is given. Assume that  does not need an input file")
    qmc.run(NULL);
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


