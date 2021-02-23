#if COMPILATION_INSTRUCTIONS
(echo "#include<"$0">" > $0x.cpp) && mpicxx -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_WINDOW $0x.cpp -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_REMOTE_HPP
#define BOOST_MPI3_REMOTE_HPP


#ifdef _TEST_BOOST_MPI3_REMOTE_HPP
int main(){
}
#endif
#endif

