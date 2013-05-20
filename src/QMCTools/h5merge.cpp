#include "Utilities/RandomGenerator.h"
#include <vector>
using namespace std;
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
  string ofile(argv[1]);
  while(ic<argc)
  {
    string w(argv[ic]);
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
  cout << "Number of processors = " << np << endl;
  //qmcplusplus::HDFWalkerMerger merger(argv[1],atoi(argv[2]));
  //qmcplusplus::HDFWalkerMerger merger(argv[1],np);
  //merger.merge();
  return 0;
}
