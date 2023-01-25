//#ifdef COMPILATION_INSTRUCTIONS
//(echo '#include"'$0'"'>$0.cpp)&&mpic++ $0 -o $0x -lfftw3_mpi -lfftw3 -lboost_unit_test_framework&&mpirun -n 4 $0x&&rm $0x $0.cpp;exit
//#endif
// Copyright 2018-2021 Alfredo A. Correa

#include<fftw3-mpi.h>

#include<complex>
#include "../main.hpp"

const ptrdiff_t N0 = 16, N1 = 16;

std::complex<double> myf_x(ptrdiff_t i, ptrdiff_t j){
	return {cos(2*M_PI*i/N1), sin(2*M_PI*i/N1)};
}

std::complex<double> myf_y(ptrdiff_t i, ptrdiff_t j){
	return {cos(2*M_PI*j/N1), sin(2*M_PI*j/N1)};
}

namespace mpi3 = boost::mpi3;

int mpi3::main(int argc, char **argv, mpi3::communicator world) {
	fftw_mpi_init();

	/* get local data size and allocate */
	ptrdiff_t local_n0, local_0_start;
	ptrdiff_t alloc_local = fftw_mpi_local_size_2d(N0, N1, &world, &local_n0, &local_0_start);
	assert( alloc_local >= local_n0*N1 and alloc_local%4==0 ); // seems to be multiple of 4 (number of complexs)
	auto data = reinterpret_cast<std::complex<double>*>(fftw_alloc_complex( alloc_local )); 

	/* create plan for in-place forward DFT */
	fftw_plan plan = fftw_mpi_plan_dft_2d(
		N0, N1, reinterpret_cast<fftw_complex*>(data), reinterpret_cast<fftw_complex*>(data), 
		&world, // pass an old fashion comm handle to legacy libraries
		FFTW_FORWARD, FFTW_ESTIMATE
	);

	/* initialize data to some function my_function(x,y) */
	for(ptrdiff_t i = 0; i != local_n0; ++i)
		for(ptrdiff_t j = 0; j != N1; ++j)
			data[i*N1 + j] = myf_x(local_0_start + i, j);

	/* compute transforms, in-place, as many times as desired */
	fftw_execute(plan);
	{
		int i = 0, j = 0;
		if(local_0_start <= i and i < local_0_start + local_n0)
			assert( abs( data[(i - local_0_start)*N1 + j] )       < 1e-12 );
	}
	{
		int i = 1, j = 0;
		if(local_0_start <= i and i < local_0_start + local_n0)
			assert( abs( data[(i - local_0_start)*N1 + j] - 256.) < 1e-12 );
	}
	{
		int i = 0, j = 15;
		if(local_0_start <= i and i < local_0_start + local_n0)
			assert( abs( data[(i - local_0_start)*N1 + j] )       < 1e-12 );
	}
	fftw_destroy_plan(plan);

	fftw_mpi_cleanup();
	return 0;
}
