//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#include <iostream>
#include "SQD/SQDFrame.h"
#include <qapplication.h>

main(int argc, char **argv) {

  QApplication a( argc, argv );

  if(argc<2) {
    std::cerr << "Usage: sqd file-name [-nox]" << std::endl;
    std::cerr << "Use -nox if you like to skip plotting the results." << std::endl;
    return 1;
  }

  bool showplot = true;
  int i=1;
  while(i<argc) {
    string arg(argv[i]);
    if(arg == "-nox") {
      showplot = false;
    }
    i++;
  }

  SQDFrame *w = new SQDFrame;
  a.setMainWidget( w );
  w->solve(argv[1]);

  w->show();
  return a.exec();
}

void 
SQDFrame::solve(const char* fname) {

  using namespace ohmmshf;
  using namespace xmlpp;
  myParser.parse_file(fname);
  docRoot = myParser.get_document()->get_root_node(); //deleted by DomParser.

  parseXMLFile(Pot,Psi,elementType,potType,gridType,docRoot);

  HFSolver = new HartreeFock(Pot,Psi,docRoot);
  HFSolver->setRoot(fname);
  HFSolver->solve(potType,gridType,Psi.size());
  
   ofstream fout("test.dat");
   for(int ig; ig<Psi.m_grid->size(); ig++) {
     for(int orb=0; orb<Psi.NumUniqueOrb; orb++){
       fout << setw(15) << Psi(orb,ig);
     }
     fout << endl;
   }
  Psi.print(elementType);
}


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

