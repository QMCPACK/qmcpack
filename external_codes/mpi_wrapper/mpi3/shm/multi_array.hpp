#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -I${HOME}/prj -fpermissive -std=c++17 -Wall `#-Wfatal-errors` -D_TEST_BOOST_MPI3_SHM_MULTI $0x.cpp -o $0x.x -lblas && time mpirun -n 4 $0x.x $@; exit 
#&& rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_SHM_MULTI_HPP
#define BOOST_MPI3_SHM_MULTI_HPP

#include "../shm/allocator.hpp"

#endif