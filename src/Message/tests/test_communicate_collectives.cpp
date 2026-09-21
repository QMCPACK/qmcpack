//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <vector>

#include "Message/CommOperators.h"

namespace qmcplusplus
{
namespace
{
using Complex = std::complex<double>;

Complex rank_value(int rank, int element)
{
  const int value = 10 * rank + element + 1;
  return {static_cast<double>(value), static_cast<double>(-value)};
}
} // namespace

TEST_CASE("communicate_collectives_complex_pointer_bcast", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;
  std::vector<Complex> values{{-1, -1}, {-1, -1}};
  if (comm->rank() == 0)
    values = {{1, -2}, {3, -4}};

  comm->bcast(values.data(), static_cast<int>(values.size()));

  REQUIRE(values[0] == Complex{1, -2});
  REQUIRE(values[1] == Complex{3, -4});
}

TEST_CASE("communicate_collectives_complex_pointer_allgather", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;
  std::vector<Complex> send{rank_value(comm->rank(), 0), rank_value(comm->rank(), 1)};
  std::vector<Complex> receive(comm->size() * send.size(), Complex{-1, -1});

  comm->allgather(send.data(), receive.data(), static_cast<int>(send.size()));

  for (int rank = 0; rank < comm->size(); ++rank)
    for (int element = 0; element < static_cast<int>(send.size()); ++element)
      REQUIRE(receive[rank * send.size() + element] == rank_value(rank, element));
}

TEST_CASE("communicate_collectives_complex_pointer_gatherv", "[message][collectives]")
{
  Communicate* comm     = OHMMS::Controller;
  const int local_count = comm->rank() + 1;
  std::vector<Complex> send(local_count);
  for (int element = 0; element < local_count; ++element)
    send[element] = rank_value(comm->rank(), element);

  std::vector<int> counts(comm->size());
  std::vector<int> displacements(comm->size());
  int total_count = 0;
  for (int rank = 0; rank < comm->size(); ++rank)
  {
    counts[rank]        = rank + 1;
    displacements[rank] = total_count;
    total_count += counts[rank];
  }
  std::vector<Complex> receive(total_count, Complex{-1, -1});

  comm->gatherv(send.data(), receive.data(), local_count, counts, displacements);

  if (comm->rank() == 0)
    for (int rank = 0; rank < comm->size(); ++rank)
      for (int element = 0; element < counts[rank]; ++element)
        REQUIRE(receive[displacements[rank] + element] == rank_value(rank, element));
}

TEST_CASE("communicate_collectives_complex_matrix_allreduce", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;
  Matrix<Complex> values(2, 2);
  for (int row = 0; row < 2; ++row)
    for (int column = 0; column < 2; ++column)
    {
      const int value     = (comm->rank() + 1) * (2 * row + column + 1);
      values(row, column) = Complex{static_cast<double>(value), static_cast<double>(-value)};
    }

  comm->allreduce(values);

  const int rank_sum = comm->size() * (comm->size() + 1) / 2;
  for (int row = 0; row < 2; ++row)
    for (int column = 0; column < 2; ++column)
    {
      const int value = rank_sum * (2 * row + column + 1);
      REQUIRE(values(row, column) == Complex{static_cast<double>(value), static_cast<double>(-value)});
    }
}

TEST_CASE("communicate_collectives_serial_allgather_count", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;
  if (comm->size() != 1)
    return;

  std::vector<int> send{1, 2, 3};
  std::vector<int> receive{0, 0, -1};
  comm->allgather(send, receive, 2);

  REQUIRE(receive[0] == 1);
  REQUIRE(receive[1] == 2);
  REQUIRE(receive[2] == -1);
}
} // namespace qmcplusplus
