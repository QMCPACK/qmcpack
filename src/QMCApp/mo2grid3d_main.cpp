//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "QMCApp/MO2Grid3D.h"

/**@file mo2grid3d_main.cpp
 *@brief a main function to map MolecularOrbitals on 3-D numerical Orbitals
 *
 * Using MO2Grid3D as the engine to transform MolecularOrbitals
 */
int main(int argc, char **argv)
{
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv,OHMMS::Controller->mycontext());
  qmcplusplus::MO2Grid3D qmc(argc,argv);
  if(argc>1)
  {
    if(qmc.parse(argv[1]))
    {
      qmc.execute();
      //qmc.saveXml();
    }
    //xmlFreeDoc(m_doc);
  }
  else
  {
    ERRORMSG("No input file is given.")
    ERRORMSG("usage: mo2grid3d input-file")
  }
  LOGMSG("Bye")
  OHMMS::Controller->finalize();
  return 0;
}


