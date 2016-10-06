//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Utilities/RandomGenerator.h"
#include <vector>
#include "QMCTools/HDFWalkerMerger.h"
#include "Utilities/OhmmsInfo.h"

int main(int argc, char **argv)
{
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv,OHMMS::Controller->rank());
  qmcplusplus::Random.init(0,1,-1);
  if(argc<2)
  {
    std::cerr << " Usage: h5merge <rootname> -n <number-of-processor> -o[utfile] <outfile> " << std::endl;
    return 1;
  }
  int ic=0;
  int np=1;
  std::string ofile(argv[1]);
  while(ic<argc)
  {
    std::string w(argv[ic]);
    if(w.find("-n")<w.size())
    {
      np=atoi(argv[++ic]);
    }
    else
      if(w.find("-o")<w.size())
      {
        ofile=argv[++ic];
      }
    ++ic;
  }
  std::cout << "Number of processors = " << np << std::endl;
  //qmcplusplus::HDFWalkerMerger merger(argv[1],atoi(argv[2]));
  //qmcplusplus::HDFWalkerMerger merger(argv[1],np);
  //merger.merge();
  return 0;
}
