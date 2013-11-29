//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file phdf.cpp
 * @brief Test code for PHDF
 */
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include <Configuration.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <algorithm>
#include <io/hdf_archive.h>

using namespace qmcplusplus;

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  Communicate* comm=OHMMS::Controller;
  OhmmsInfo Welcome(argc,argv,comm->rank());
  Random.init(0,1,11);

  const int n=64;
  TinyVector<int,3> gcounts(n,4,3);
  TinyVector<int,3> counts(n,4,3);
  TinyVector<int,3> offsets(0,0,0);

  counts[0]=n/comm->size();
  offsets[0]=(n/comm->size())*comm->rank();

  cout << comm->rank() << " " << counts << " " << offsets << endl;

  typedef vector<double> buffer_t;
  buffer_t mydata(counts[0],comm->rank());

  hdf_archive hout(comm,true);
  hout.create("test.h5");
  hout.write(n,"size");

  //hyperslab_proxy<buffer_t,3> slab(mydata,gcounts,counts,offsets);
  //hout.write(slab,"hello");
  OHMMS::Controller->finalize();
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $
 ***************************************************************************/
