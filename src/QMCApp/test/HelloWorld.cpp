#include "Configuration.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/SimpleParser.h"
#include "OhmmsData/FileUtility.h"
#include "Platforms/sysutil.h"
#include "Platforms/devices.h"
#include "OhmmsApp/ProjectData.h"
#include "QMCApp/QMCMain.h"
#include "qmc_common.h"


int main(int argc, char **argv)
{
  using namespace qmcplusplus;
  OHMMS::Controller->initialize(argc,argv);


  OHMMS::Controller->finalize();
  return 0;
}


