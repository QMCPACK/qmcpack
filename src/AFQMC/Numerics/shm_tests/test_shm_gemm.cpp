//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#undef NDEBUG

//#include "catch.hpp"
#include "Configuration.h"

#include <vector>
#include <iostream>

#include "mpi3/communicator.hpp"

#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "Utilities/FairDivide.h"
#include "AFQMC/Utilities/myTimer.h"

using boost::multi::array;
using boost::multi::array_ref;
using std::cout;
using std::endl;
using std::vector;

using namespace qmcplusplus;
using namespace afqmc;

void timing_shm_blas(int c)
{
  using Type         = double;
  using communicator = boost::mpi3::communicator;
  using shm_Alloc    = shared_allocator<Type>;
  using Matrix       = boost::multi::array_ref<Type, 2>;

  myTimer Timer;
  int n0 = 64, npower = 6, nmax = n0 * std::pow(2, npower - 1);
  int ntimes = 5;

  auto world = boost::mpi3::environment::get_world_instance();
  auto node  = world.split_shared();

  int memory_needs = nmax * (c * c * nmax + 2 * c * nmax);
  boost::multi::array<Type, 1, shared_allocator<Type>> buff(iextensions<1u>{memory_needs},
                                                            shared_allocator<Type>{node});

  std::vector<std::pair<int, int>> pairs;
  {
    int imax = std::ceil(std::sqrt(node.size())) + 1;
    pairs.push_back({node.size(), 1});
    for (int i = 2; i <= imax; i++)
      if (node.size() % i == 0)
        pairs.push_back({node.size() / i, i});
  }

  for (int p = 0; p < npower; p++)
  {
    int n  = n0 * std::pow(2, p);
    int nu = c * n;
    Matrix A(buff.data(), {nu, nu});
    Matrix B(buff.data() + nu * nu, {nu, n});
    Matrix C(buff.data() + nu * nu + nu * n, {nu, n});

    for (auto& v : pairs)
    {
      int nr = v.first;
      int nc = v.second;

      int myrow = node.rank() / nc;
      int mycol = node.rank() % nc;

      // split over rows
      int r0, rN, c0, cN;
      std::tie(r0, rN) = FairDivideBoundary(myrow, nu, nr);
      std::tie(c0, cN) = FairDivideBoundary(mycol, n, nc);


      ma::product(A.sliced(r0, rN), B(B.extension(0), {c0, cN}), C({r0, rN}, {c0, cN}));

      Timer.reset("Gen");
      node.barrier();
      Timer.start("Gen");
      for (int t = 0; t < ntimes; t++)
      {
        ma::product(A.sliced(r0, rN), B(B.extension(0), {c0, cN}), C({r0, rN}, {c0, cN}));
        node.barrier();
      }
      Timer.stop("Gen");
      if (node.root())
        cout << " [" << nr << "][" << nc << "]: " << nu << " " << n << " " << Timer.total("Gen") / ntimes << endl;
    }
  }
  if (node.root())
    cout << endl;
}

int main(int argc, char* argv[])
{
  boost::mpi3::communicator world{MPI_COMM_WORLD};
  auto world = boost::mpi3::environment::get_world_instance();
  int c      = 10;
  if (argc > 1)
    c = atoi(argv[1]);
  timing_shm_blas(c);
  OHMMS::Controller->finalize();
}
