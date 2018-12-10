#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -Wall -Wextra -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_EQUALITY $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x ddd $@ && rm -f $0x.x $0x.cpp; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#ifndef BOOST_MPI3_EQUALITY_HPP
#define BOOST_MPI3_EQUALITY_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

namespace boost{
namespace mpi3{

enum equality{
	identical = MPI_IDENT, 
	congruent = MPI_CONGRUENT, 
	similar = MPI_SIMILAR, 
	unequal = MPI_UNEQUAL
};

}}

#ifdef _TEST_BOOST_MPI3_EQUALITY

#include<iostream>
#include<chrono>

namespace mpi3 = boost::mpi3;
using std::cout;

int main(int, char*[]){
	mpi3::equality q = mpi3::identical;
	cout << q << std::endl;
	q = mpi3::congruent;
	cout << q << std::endl;
	q = mpi3::similar;
	cout << q << std::endl;
	q = mpi3::unequal;
	cout << q << std::endl;
}

#endif
#endif

