//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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

  std::cout << comm->rank() << " " << counts << " " << offsets << std::endl;

  typedef std::vector<double> buffer_t;
  buffer_t mydata(counts[0],comm->rank());

  hdf_archive hout(comm,true);
  hout.create("test.h5");
  hout.write(n,"size");

  //hyperslab_proxy<buffer_t,3> slab(mydata,gcounts,counts,offsets);
  //hout.write(slab,"hello");
  OHMMS::Controller->finalize();
  return 0;
}

