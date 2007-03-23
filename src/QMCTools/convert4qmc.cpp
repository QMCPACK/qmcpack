#include "QMCTools/CasinoParser.h"
#include "QMCTools/GaussianFCHKParser.h"
#include "QMCTools/GamesXmlParser.h"
#include "QMCTools/BParser.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "OhmmsData/FileUtility.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"

int main(int argc, char **argv) {

  if(argc<2) {
    std::cout << "Usage: convert [-gaussian|-casino|-gamesxml] filename ";
    std::cout << "[-hdf5 -psi_tag psi0 -ion_tag ion0 -gridtype log|log0|linear -first ri -last rf -size npts]" << std::endl;
    std::cout << "Defaults : -gridtype log -first 1e-6 -last 100 -size 1001" << std::endl;
    std::cout << "When the input format is missing, the  extension of filename is used to determine the parser " << std::endl;
    std::cout << " *.Fchk -> gaussian; *.data -> casino; *.xml -> gamesxml" << std::endl;
    return 1;
  }


  OHMMS::Controller->initialize(argc,argv);

  OhmmsInfo welcome(argc,argv,OHMMS::Controller->mycontext());
  Random.init(0,1,-1);

  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right,std::ios::adjustfield);
  std::cout.precision(12);

  QMCGaussianParserBase::init();

  QMCGaussianParserBase *parser=0;
  int iargc=0;
  string in_file(argv[1]);
  string psi_tag("psi0");
  string ion_tag("ion0");
  bool usehdf5=false;
  while(iargc<argc) {
    std::string a(argv[iargc]);
    if(a == "-gaussian") {
      parser = new GaussianFCHKParser(argc,argv);
      in_file =argv[++iargc];
    }
    else if(a == "-gamesxml") {
      parser = new GamesXmlParser(argc,argv);
      in_file =argv[++iargc];
    }
    else if(a == "-casino") {
      parser = new CasinoParser(argc,argv);
      in_file =argv[++iargc];
    }
    else if(a == "-b") {
      parser = new BParser(argc,argv);
      in_file =argv[++iargc];
    } else if(a == "-hdf5") {
      usehdf5=true;
    } else if(a == "-psi_tag") {
      psi_tag=argv[++iargc];
    } else if(a == "-ion_tag") {
      ion_tag=argv[++iargc];
    }
    ++iargc;
  }

  //Failed to create a parser. Try with the extension
  if(parser == 0) {
    string ext= getExtension(in_file);
    if(ext == "data") {
      WARNMSG("Creating CasinoParser")
      parser = new CasinoParser(argc,argv);
    } else if(ext == "Fchk") {
      WARNMSG("Creating GaussianFCHKParser")
      parser = new GaussianFCHKParser(argc,argv);
    } else if(ext == "xml") {
      WARNMSG("Creating GamesXmlParser")
      parser = new GamesXmlParser(argc,argv);
    } else if(ext == "10") {
      WARNMSG("Creating BParser")
      parser = new BParser(argc,argv);
    }
  }

  parser->UseHDF5=usehdf5;
  parser->IonSystem.setName(ion_tag);
  parser->parse(in_file);
  parser->dump(psi_tag, ion_tag);
  return 0;
}

/*
int main(int argc, char **argv) {

  char buffer[200];
  std::string _txt;
  std::string _data;
  double _t;
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right,std::ios::adjustfield);
  std::cout.precision(12);
  while(std::cin.getline( buffer, sizeof ( buffer ) ) ){
    std::istringstream stream(buffer);
    if(isdigit(buffer[1])) {
      while(stream>>_t) {
        std::cout << std::setw(21) << _t ;
      }
      std::cout << std::endl;
    } else {
      if(stream>>_t) {
        std::cout << "probably numbers " << buffer << std::endl;
      } else {
        std::cout << "probably statement " << buffer << std::endl;
      }
    }
  }
}
*/
