//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "QMCTools/CasinoParser.h"
#include "QMCTools/GaussianFCHKParser.h"
#include "QMCTools/GamesXmlParser.h"
#include "QMCTools/GamesAsciiParser.h"
#include "QMCTools/VSVBParser.h"
#include "QMCTools/QPParser.h"
#include "QMCTools/GamesFMOParser.h"
#include "QMCTools/BParser.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "OhmmsData/FileUtility.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"

int main(int argc, char **argv)
{
  if(argc<2)
  {
    std::cout << "Usage: convert [-gaussian|-casino|-gamesxml|-gamessAscii|-gamessFMO |-VSVB| -QP] filename ";
     std::cout << "[-nojastrow -hdf5 -psi_tag psi0 -ion_tag ion0 -gridtype log|log0|linear -first ri -last rf -size npts -ci file.out -threshold cimin -TargetState state_number -NaturalOrbitals NumToRead -add3BodyJ -prefix title -addCusp]"
              << std::endl;
    std::cout << "Defaults : -gridtype log -first 1e-6 -last 100 -size 1001 -ci required -threshold 0.01 -TargetState 0 -prefix sample" << std::endl;
    std::cout << "When the input format is missing, the  extension of filename is used to determine the parser " << std::endl;
    std::cout << " *.Fchk -> gaussian; *.out -> gamessAscii; *.data -> casino; *.xml -> gamesxml" << std::endl;
    return 1;
  }
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv,OHMMS::Controller->rank());
  Random.init(0,1,-1);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right,std::ios::adjustfield);
  std::cout.precision(12);
  QMCGaussianParserBase::init();
  QMCGaussianParserBase *parser=0;
  int iargc=0;
  std::string in_file(argv[1]);
  std::string punch_file;
  std::string psi_tag("psi0");
  std::string ion_tag("ion0");
  std::string prefix("sample");
  int TargetState=0;
  bool usehdf5=false;
  bool ci=false,zeroCI=false,orderByExcitation=false,VSVB=false, fmo=false,addCusp=false;
  double thres=0.01;
  int readNO=0; // if > 0, read Natural Orbitals from gamess output
  int readGuess=0; // if > 0, read Initial Guess from gamess output
  while(iargc<argc)
  {
    std::string a(argv[iargc]);
    if(a == "-gaussian")
    {
      parser = new GaussianFCHKParser(argc,argv);
      in_file =argv[++iargc];
    }
    else if(a == "-gamesxml")
    {
      parser = new GamesXmlParser(argc,argv);
      in_file =argv[++iargc];
    }
    else if(a == "-gamessAscii")
    {
      parser = new GamesAsciiParser(argc,argv);
      in_file =argv[++iargc];
    }
    else if(a == "-QP")
    {
      parser = new QPParser(argc,argv);
      in_file =argv[++iargc];
    }
    else if(a == "-VSVB")
    {
      parser = new VSVBParser(argc,argv);
      in_file =argv[++iargc];
      VSVB=true;
    }
    else if(a == "-gamessFMO")
    {
      parser = new GamesFMOParser(argc,argv);
      in_file =argv[++iargc];
      fmo=true;
    }
    else if(a == "-casino")
    {
      parser = new CasinoParser(argc,argv);
      in_file =argv[++iargc];
    }
    else if(a == "-b")
    {
      parser = new BParser(argc,argv);
      in_file =argv[++iargc];
    }
    else if(a == "-hdf5")
    {
      usehdf5=true;
    }
    else if(a == "-psi_tag")
    {
      psi_tag=argv[++iargc];
    }
    else if(a == "-ion_tag")
    {
      ion_tag=argv[++iargc];
    }
    else if(a == "-prefix")
    {
      prefix=argv[++iargc];
    }
    else if(a == "-ci")
    {
      ci=true;
      punch_file = argv[++iargc];
    }
    else if(a == "-addCusp" )
    {
      addCusp = true;
    }
    else if(a == "-threshold" )
    {
      thres = atof(argv[++iargc]);
    }
    else if(a == "-TargetState" )
    {
      TargetState = atoi(argv[++iargc]);
    }
    else if(a == "-NaturalOrbitals")
    {
      readNO = atoi(argv[++iargc]);
    }
    else if(a == "-readInitialGuess")
    {
      readGuess = atoi(argv[++iargc]);
    }
    else if(a == "-zeroCI")
    {
      zeroCI = true;
    }
    else if(a == "-orderByExcitation")
    {
      orderByExcitation = true;
    }
    else if(a == "-cutoff")
    {
      orderByExcitation = true;
    }
    ++iargc;
  }
  if(readNO > 0 && readGuess > 0)
  {
    std::cerr <<"Can only use one of: -NaturalOrbitals or -readInitialGuess. \n";
    abort();
  }
  //Failed to create a parser. Try with the extension
  if(parser == 0)
  {
    std::string ext= getExtension(in_file);
    if(ext == "data")
    {
      WARNMSG("Creating CasinoParser")
      parser = new CasinoParser(argc,argv);
    }
    else if(ext == "Fchk")
    {
      WARNMSG("Creating GaussianFCHKParser")
      parser = new GaussianFCHKParser(argc,argv);
    }
    else if(ext == "xml")
    {
      WARNMSG("Creating GamesXmlParser")
      parser = new GamesXmlParser(argc,argv);
    }
    else if(ext == "10")
    {
      WARNMSG("Creating BParser")
      parser = new BParser(argc,argv);
    }
    else if(ext == "out")
    {
      WARNMSG("Creating GamesAsciiParser")
      parser = new GamesAsciiParser(argc,argv);
    }
    else
    {
      std::cerr << "Unknown extension: " << ext << std::endl;
      exit(1);
    }
  }
  if(fmo)
  {
    parser->Title=prefix;
    parser->UseHDF5=usehdf5;
    parser->DoCusp=addCusp;
    parser->parse(in_file);

  }
  else
  {
    parser->Title=prefix;
    parser->DoCusp=addCusp;
    parser->UseHDF5=usehdf5;
    if (usehdf5)
      parser->h5file=parser->Title+".eig.h5";
    parser->IonSystem.setName(ion_tag);
    parser->multideterminant=ci;
    parser->ci_threshold=thres;
    parser->target_state=TargetState;
    parser->readNO=readNO;
    parser->orderByExcitation=orderByExcitation;
    parser->zeroCI = zeroCI;
    parser->readGuess=readGuess;
    parser->outputFile=punch_file;
    parser->VSVB=VSVB;
    parser->parse(in_file);
    parser->dump(psi_tag, ion_tag);
    OHMMS::Controller->finalize();
    

  }
  return 0;
}

