#if COMPILATION_INSTRUCTIONS
mpic++ -I$HOME/prj/alf $0 -o $0x -lfftw3 -lfftw3_mpi&&time mpirun -n 4 $0x&&rm $0x;exit
#endif

#include "../../../fftw/mpi.hpp"

#include<boost/mpi3/main.hpp>
#include<boost/mpi3/environment.hpp>
#include<boost/mpi3/ostream.hpp>
#include "../../../fftw.hpp"

namespace mpi3  = boost::mpi3;
namespace multi = boost::multi;

int mpi3::main(int, char*[], mpi3::communicator world){
	multi::fftw::mpi::environment fenv;

	multi::fftw::mpi::array<std::complex<double>, 2> G({41, 321}, world);

	if(auto x = G.local_cutout().extensions())
		for(auto i : std::get<0>(x))
			for(auto j : std::get<1>(x))
				G.local_cutout()[i][j] = std::complex<double>(i + j, i + 2*j);
	
	multi::array<std::complex<double>, 2> L = G; // world replicas
	assert( L == G );
	
	using multi::fftw::dft_forward;

	dft_forward(L, L); // dft in replicas
	dft_forward(G, G);

	if(auto x = G.local_cutout().extensions())
		for(auto i : std::get<0>(x))
			for(auto j : std::get<1>(x))
				if(not(std::abs(G.local_cutout()[i][j] - L[i][j]) < 1e-8)) std::cout<< std::abs(G.local_cutout()[i][j] - L[i][j]) << std::endl;

	return 0;
}

