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
    
    



#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "Configuration.h"
#include "Utilities/RandomGenerator.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Message/CommunicateGroup.h"
using namespace qmcplusplus;

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv,OHMMS::Controller->rank());
  Random.init(OHMMS::Controller->rank(),OHMMS::Controller->size(),-1);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right,std::ios::adjustfield);
  std::cout.precision(12);
  double sumL=0.0,sumG=0.0;
  for(int i=0; i<1000; i++)
    sumL += Random();
  int ndiv=atoi(argv[1]);
  std::cout << "Group = " << ndiv << std::endl;
  Communicate newComm(*OHMMS::Controller,ndiv);
  std::cout << OHMMS::Controller->rank() << " Rank = "
            << newComm.rank() << " Size = " << newComm.size() << " " << sumL << std::endl;
  sumG=sumL;
  newComm.allreduce(sumG);
  //MPI_Allreduce(&(sumL), &(sumG), 1, MPI_DOUBLE, MPI_SUM, newComm.getMPI());
  //int p=OHMMS::Controller->mycontext()/ndiv;
  //int q=OHMMS::Controller->mycontext()%ndiv;
  //MPI_Comm row;
  //MPI_Comm_split(OHMMS::Controller->getID(),p,q,&row);
  //int new_rank, new_size;
  //MPI_Comm_rank(row,&new_rank);
  //MPI_Comm_size(row,&new_size);
  //std::cout << OHMMS::Controller->mycontext() << " Rank = " << new_rank << " Size = " << new_size << " " << sumL << std::endl;
  //MPI_Allreduce(&(sumL), &(sumG), 1, MPI_DOUBLE, MPI_SUM, row);
  std::cout << OHMMS::Controller->rank() << " Local sum = " << sumL << " Global sum " << sumG << std::endl;
  //const int ndims=2;
  //int dims[]={2,2};
  //int coords[ndims];
  //bool periods[]={false,false};
  //MPI_Comm what=c.Get_mpi();
  //OOMPI_Cart_comm cart(OOMPI_COMM_WORLD,ndims,dims,periods);
  //cart.Coords(coords);
  //std::cout << " Rank = " << cart.Rank() << " " << coords[0] << " " << coords[1] << std::endl;
  //
  //std::cout << "World communicator " << OHMMS::Controller->getID() <<std::endl;
  //OOMPI_Intra_comm c(OHMMS::Controller->getID());
  //THIS DOES NOT WORK NO IDEA
  //c.Split(p,q);
  //std::cout << c.Get_mpi() << " Rank = " << c.Rank() << " Size = " << c.Size() << std::endl;
  //alternatively, use group
  //CommunicateGroup testGroup(*(OHMMS::Controller),ndiv);
  //sumL=0.0;
  //for(int i=0; i<1000; i++) sumL += Random();
  //std::cout << OHMMS::Controller->mycontext() << " Rank = " << testGroup.mycontext() << " Size = " << testGroup.ncontexts() << " " << sumL << std::endl;
  //MPI_Allreduce(&(sumL), &(sumG), 1, MPI_DOUBLE, MPI_SUM, testGroup.getID());
  //std::cout << OHMMS::Controller->mycontext() << " Local sum = " << sumL << " Global sum " << sumG << std::endl;
  OHMMS::Controller->finalize();
  return 0;
}
