//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <iostream>
#include "SQD/SQDFrame.h"
#include <qapplication.h>

main(int argc, char **argv)
{
  QApplication a( argc, argv );
  if(argc<2)
  {
    std::cerr << "Usage: sqd file-name [-nox]" << std::endl;
    std::cerr << "Use -nox if you like to skip plotting the results." << std::endl;
    return 1;
  }
  bool showplot = true;
  int i=1;
  while(i<argc)
  {
    std::string arg(argv[i]);
    if(arg == "-nox")
    {
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
SQDFrame::solve(const char* fname)
{
  using namespace ohmmshf;
  using namespace xmlpp;
  myParser.parse_file(fname);
  docRoot = myParser.get_document()->get_root_node(); //deleted by DomParser.
  parseXMLFile(Pot,Psi,elementType,potType,gridType,docRoot);
  HFSolver = new HartreeFock(Pot,Psi,docRoot);
  HFSolver->setRoot(fname);
  HFSolver->solve(potType,gridType,Psi.size());
  std::ofstream fout("test.dat");
  for(int ig; ig<Psi.m_grid->size(); ig++)
  {
    for(int orb=0; orb<Psi.NumUniqueOrb; orb++)
    {
      fout << std::setw(15) << Psi(orb,ig);
    }
    fout << std::endl;
  }
  Psi.print(elementType);
}



