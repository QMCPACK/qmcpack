#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
using namespace qmcplusplus;

int main(int argc, char** argv) {

  OHMMS::Controller->initialize(argc,argv);

  OhmmsInfo welcome(argc,argv,OHMMS::Controller->mycontext());
  Random.init(OHMMS::Controller->mycontext(),OHMMS::Controller->ncontexts(),-1);

  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right,std::ios::adjustfield);
  std::cout.precision(12);

  const int ndims=2;
  int dims[]={2,2};
  int coords[ndims];
  bool periods[]={false,false};

  int ndiv=atoi(argv[1]);

  int p=OHMMS::Controller->mycontext()/ndiv;
  int q=OHMMS::Controller->mycontext()%ndiv;

  //MPI_Comm what=c.Get_mpi();
  //OOMPI_Cart_comm cart(OOMPI_COMM_WORLD,ndims,dims,periods);
  //cart.Coords(coords);

  //std::cout << " Rank = " << cart.Rank() << " " << coords[0] << " " << coords[1] << std::endl;
  //
  //std::cout << "World communicator " << OHMMS::Controller->getID() <<std::endl;
  //OOMPI_Intra_comm c(OHMMS::Controller->getID());

  //c.Split(p,q);
  //std::cout << c.Get_mpi() << " Rank = " << c.Rank() << " Size = " << c.Size() << std::endl;

  MPI_Comm row;
  MPI_Comm_split(OHMMS::Controller->getID(),p,q,&row);
  int new_rank, new_size;
  MPI_Comm_rank(row,&new_rank);
  MPI_Comm_size(row,&new_size);

  double sumL=0.0,sumG=0.0;
  for(int i=0; i<1000; i++) sumL += Random();

  std::cout << OHMMS::Controller->mycontext() << " Rank = " << new_rank << " Size = " << new_size << " " << sumL << std::endl;

  MPI_Allreduce(&(sumL), &(sumG), 1, MPI_DOUBLE, MPI_SUM, row);

  std::cout << OHMMS::Controller->mycontext() << " Local sum = " << sumL << " Global sum " << sumG << std::endl;

  OHMMS::Controller->finalize();
}
