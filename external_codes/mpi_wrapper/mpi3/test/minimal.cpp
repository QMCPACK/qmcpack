#if COMPILATION_INSTRUCTIONS
$CXX `mpic++ -showme:compile|sed 's/-pthread/ /g'` $0 -o $0x `mpic++ -showme:link|sed 's/-pthread/ /g'`&&
#(mpirun -n 3 valgrind --leak-check=full --show-reachable=yes --error-limit=no                                 --suppressions=$0.openmpi.supp $0x)&&rm $0x;exit
mpic++    -O3 -Wall -Wextra $0 -o $0x&&(mpirun -n 3 valgrind --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all $0x 2>&1|grep -v '==' > $0.openmpi.supp    )&&
mpic++ -g -O3 -Wall -Wextra $0 -o $0x&&(mpirun -n 3 valgrind --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all $0x 2>&1|grep -v '==' > $0-g.openmpi.supp  )&&
mpic++ -g     -Wall -Wextra $0 -o $0x&&(mpirun -n 3 valgrind --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all $0x 2>&1|grep -v '==' > $0-O0.openmpi.supp  )&&
cat $0-g.openmpi.supp $0-O0.openmpi.supp >> $0.openmpi.supp&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020
#include "../../mpi3/environment.hpp"

namespace mpi3 = boost::mpi3;

int main(int argc, char* argv[]){
//	mpi3::environment env(argc, argv);
//	MPI_Init(&argc, &argv);
//	MPI_Finalize();

	MPI_Init(&argc, &argv);
	MPI_Comm dup_comm_world;
	MPI_Comm_dup( MPI_COMM_WORLD, &dup_comm_world );
	MPI_Comm_free(&dup_comm_world);
	MPI_Finalize();

}

