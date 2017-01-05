///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/random/random.cpp
///
/// \brief   implementation for functions and objects associated with random numbers
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <formic/utils/random.h>
#include <formic/utils/mpi_interface.h>

// initialize the global linear congruential random number generator
formic::LinearCongruentialGenerator formic::global_lcg;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   sets the random seed on each process
///
/// \param[in]     seed     the random seed
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::set_seed(unsigned int seed) {

  // get MPI info
  const int nproc  = formic::mpi::size();
  const int myrank = formic::mpi::rank();

  // temporarily initialize using the provided seed
  formic::global_lcg.set_seed(seed);

  // produce a different seed for each process
  std::vector<unsigned int> seeds(nproc);
  for (int i = 0; i < nproc; i++) {
    unsigned int temp;
    for (int j = 0; j < 1000000; j++) temp = formic::global_lcg(); // separate the seeds by a million random numbers
    seeds.at(i) = formic::global_lcg();
  }

  // give every process a different seed
  formic::mpi::scatter(&seeds.at(0), &seed, 1, 0);

  // initialize using this process's seed
  formic::global_lcg.set_seed(seed);

}
