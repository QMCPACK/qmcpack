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
    
    



#include <string>
#include <sstream>
#include <mpi.h>

#include <map>
#include <iostream>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>

template<class T>
void mpi_send_object(const T& obj, int proc, int tag, const MPI_Comm& comm)
{
  std::ostringstream oss;
  boost::archive::binary_oarchive oa(oss);
  oa << obj;
  std::string str(oss.str());
  int size = str.size();
  MPI_Ssend(&size, 1, MPI_INT, proc, tag, comm);
  MPI_Ssend(&str[0], size, MPI_CHAR, proc, tag, comm);
}

template<class T>
void mpi_recv_object(T& obj, int proc, int tag, const MPI_Comm& comm)
{
  int size;
  MPI_Status status;
  MPI_Recv(&size, 1, MPI_INT, proc, tag, comm, &status);
  std::string buf(size, ' ');
  MPI_Recv(&buf[0], size, MPI_CHAR, proc, tag, comm, &status);
  std::istringstream iss(buf);
  boost::archive::binary_iarchive ia(iss);
  ia >> obj;
}


int main(int nargs, char** argv)
{
  MPI_Init(&nargs, &argv);
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  if (comm_size!=2)
  {
    std::cerr << "communicator size is " << comm_size
              << "; it must be 2 " << std::endl;
    std::exit(2);
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank==0)
  {
    std::map<int, std::string> m;
    m[1] = "h1";
    m[10] = "h10";
    m[56] = "56";
    m[500] = "560s";
    mpi_send_object(m, 1, 0, MPI_COMM_WORLD);
  }
  else
  {
    std::map<int, std::string> m;
    mpi_recv_object(m, 0, 0, MPI_COMM_WORLD);
    for (std::map<int, std::string>::iterator it=m.begin(); it!=m.end(); ++it)
    {
      std::cout << "m[" << it->first << "]=" << it->second << std::endl;
    }
  }
  MPI_Finalize();
  return 0;
}
