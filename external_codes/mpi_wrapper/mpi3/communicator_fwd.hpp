#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -g -std=c++14 `#-Wfatal-errors` -D_TEST_BOOST_MPI3_COMMUNICATOR_FWD $0x.cpp -o $0x.x && time mpirun -np 1 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_COMMUNICATOR_FWD_HPP
#define BOOST_MPI3_COMMUNICATOR_FWD_HPP

namespace boost{
namespace mpi3{

struct communicator;

}}

#ifdef _TEST_BOOST_MPI3_COMMUNICATOR_FWD

#include<iostream>
using std::cout;

int main(){}

#endif
#endif


