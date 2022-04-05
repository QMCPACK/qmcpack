//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifdef __linux__
#include <sched.h>
#endif
#include <iostream>
#include <sstream>
#include "config.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include "Concurrency/OpenMP.h"

namespace
{
/*=======================================*/
/* routine to return the core ID         */
/*=======================================*/
int get_core()
{
#ifdef __linux__
  int cpuid = sched_getcpu();
  return cpuid;
#else
  return -1;
#endif
}

/*=======================================*/
/* routine to return the HW thread ID    */
/*=======================================*/
int get_hwthread() { return -1; }

} // namespace

int main()
{
  int rank = 0, world_size = 1;
#ifdef HAVE_MPI
  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif

  bool hwthread_id_supported = (get_hwthread() != -1);

  if (rank == 0)
    std::cout << "OpenMP OMP_MAX_ACTIVE_LEVELS = " << omp_get_max_active_levels() << std::endl;

  for (int l_rank = 0; l_rank < world_size; l_rank++)
  {
    if (l_rank == rank)
    {
#pragma omp parallel
      {
        int L1_thread_id  = omp_get_thread_num();
        int L1_num_thread = omp_get_num_threads();
#pragma omp parallel
        {
#pragma omp master
          if (L1_thread_id == 0 and l_rank == 0)
            std::cout << "OpenMP enables " << L1_num_thread << " 1st level threads, "
                      << "and " << omp_get_num_threads() << " 2nd level threads." << std::endl
                      << std::endl;
        }
#pragma omp barrier

        for (int i = 0; i < L1_num_thread; i++)
        {
          if (L1_thread_id == i)
          {
#pragma omp parallel
            {
              int L2_thread_id  = omp_get_thread_num();
              int L2_num_thread = omp_get_num_threads();
              for (int j = 0; j < L2_num_thread; j++)
              {
                if (L2_thread_id == j)
                {
                  std::ostringstream o;
                  o << "MPI rank " << rank << " L1 tid " << L1_thread_id << " L2 tid " << L2_thread_id
                    << " is placed on Core ID " << get_core();
                  if (hwthread_id_supported)
                    o << " HW thread ID " << get_hwthread();
                  o << std::endl;
                  std::cout << o.str();
                }
#pragma omp barrier
              }
            }
          }
#pragma omp barrier
        }
      }
    }
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
