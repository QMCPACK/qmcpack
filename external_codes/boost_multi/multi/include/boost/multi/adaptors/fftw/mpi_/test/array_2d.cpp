// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

//#include "../../../fftw/mpi.hpp"

#include<mpi3/main.hpp>
#include<mpi3/environment.hpp>
#include<mpi3/ostream.hpp>
#include "../../../fftw.hpp"

namespace mpi3  = boost::mpi3;
namespace multi = boost::multi;

int mpi3::main(int, char**, mpi3::communicator /*world*/){
//	multi::fftw::mpi::environment fenv;

//  multi::fftw::mpi::array<std::complex<double>, 2> G({41, 321}, world);

#if 0
	if(auto x = G.local_cutout().extensions())
		for(auto i : std::get<0>(x))
			for(auto j : std::get<1>(x))
				G.local_cutout()[i][j] = std::complex<double>(i + j, i + 2*j);

	multi::array<std::complex<double>, 2> L = G;  // world replicas
	assert( L == G );

	using multi::fftw::dft_forward;

	dft_forward(L, L);  // dft in replicas
	dft_forward(G, G);

	if(auto x = G.local_cutout().extensions())
		for(auto i : std::get<0>(x))
			for(auto j : std::get<1>(x))
				if(not(std::abs(G.local_cutout()[i][j] - L[i][j]) < 1e-8)) std::cout<< std::abs(G.local_cutout()[i][j] - L[i][j]) << std::endl;
#endif
	return 0;
}

