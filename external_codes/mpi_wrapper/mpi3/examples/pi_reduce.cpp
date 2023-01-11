#if COMPILATION_INSTRUCTIONS
OMPI_CXX=$CXXX mpicxx $CXXFLAGS $0 -o $0x&&mpirun --oversubscribe -n 6 $0x&&rm $0x;exit
#endif

#include "../examples/../process.hpp"
#include "../examples/../main_environment.hpp"

#include <random>

namespace mpi3 = boost::mpi3;

// this program is based on https://computing.llnl.gov/tutorials/mpi/samples/C/mpi_pi_reduce.c



double pi_approx(int n_samples){

	std::random_device rd;
	static std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-1.0, +1.0);

	int score = 0;
	for(int i = 0; i != n_samples; ++i){
		double x = dis(gen);
		double y = dis(gen);
		if(x*x + y*y <= 1.0) ++score;
	}
	return 4.*score/static_cast<double>(n_samples); // pi
} 

int mpi3::main(int, char*[], mpi3::environment& menv){
	auto world = menv.world();

	constexpr int n_samples = 50000000;

	for(int s = 0; s != world.size(); ++s){
		auto t0 = mpi3::wall_time();
		if(auto continent = (world <= s)){

			double pi = (continent += pi_approx(n_samples/continent.size())/continent.size()); 

			if(world.root()){
				std::cout
					<<"After "<< n_samples <<" throws "
					<<"(in "<< continent.size() <<" procs), "
					<<"average value of pi = "   << pi	
					<<" (vs. real value of pi = "<< M_PI <<")"
					<<" (time "<< (mpi3::wall_time() - t0).count() <<" sec)"
					<<std::endl
				;
			}
		}
	}
	return 0;
}

