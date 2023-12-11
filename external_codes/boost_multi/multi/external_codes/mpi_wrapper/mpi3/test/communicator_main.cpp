// #if COMPILATION_INSTRUCTIONS
// mpic++    -O3 -Wall -Wextra $0 -o $0x&&(mpirun -n 3 valgrind --leak-check=full --show-reachable=yes --error-limit=no                                 --suppressions=$0.openmpi.supp $0x)&&rm $0x;exit
// #mpic++    -O3 -Wall -Wextra $0 -o $0x&&(mpirun -n 3 valgrind --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all $0x 2>&1|grep -v '==' > $0.openmpi.supp    )&&
// #mpic++ -g -O3 -Wall -Wextra $0 -o $0x&&(mpirun -n 3 valgrind --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all $0x 2>&1|grep -v '==' > $0-g.openmpi.supp  )&&
// #mpic++ -g     -Wall -Wextra $0 -o $0x&&(mpirun -n 3 valgrind --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all $0x 2>&1|grep -v '==' > $0-O0.openmpi.supp  )&&
// #cat $0-g.openmpi.supp $0-O0.openmpi.supp >> $0.openmpi.supp&&rm $0x;exit
// #endif
// Copyright 2018-2022 Alfredo A. Correa

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;

auto mpi3::main(int /*argv*/, char** /*argc*/, mpi3::communicator world) -> int try {
	assert( world.size() == 3 );
	return 0;
} catch(...) {return 1;}
