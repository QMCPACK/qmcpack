#include "Utilities/RandomGenerator.h"
#include <vector>
using namespace std;
#include "QMCTools/HDFWalkerMerger.h"
#include "Utilities/OhmmsInfo.h"

int main(int argc, char **argv) {

  OHMMS::Controller->initialize(argc,argv);

  OhmmsInfo welcome(argc,argv,OHMMS::Controller->mycontext());
  APPNAMESPACE::Random.init(0,1,-1);

  if(argc<2) {
    std::cerr << " Usage: h5merge rootname number-of-processor " << std::endl;
    return 1;
  }

  qmcplusplus::HDFWalkerMerger merger(argv[1],atoi(argv[2]));
  merger.merge();
  return 0;
}
