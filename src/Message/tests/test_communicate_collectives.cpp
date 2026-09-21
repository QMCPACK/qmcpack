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

  CHECK(values[0] == Complex{1, -2});
  CHECK(values[1] == Complex{3, -4});
}

TEST_CASE("communicate_collectives_complex_pointer_allgather", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;
  std::vector<Complex> send{rank_value(comm->rank(), 0), rank_value(comm->rank(), 1)};
  std::vector<Complex> receive(comm->size() * send.size(), Complex{-1, -1});

  comm->allgather(send.data(), receive.data(), static_cast<int>(send.size()));

  for (int rank = 0; rank < comm->size(); ++rank)
    for (int element = 0; element < static_cast<int>(send.size()); ++element)
      CHECK(receive[rank * send.size() + element] == rank_value(rank, element));
}

TEST_CASE("communicate_collectives_complex_pointer_gatherv", "[message][collectives]")
{
  Communicate* comm     = OHMMS::Controller;
  const int local_count = comm->rank() + 1;
  std::vector<double> send(local_count);
  for (int element = 0; element < local_count; ++element)
    send[element] = 10 * comm->rank() + element + 1;

  std::vector<int> counts(comm->size());
  std::vector<int> displacements(comm->size());
  int total_count = 0;
  for (int rank = 0; rank < comm->size(); ++rank)
  {
    counts[rank]        = rank + 1;
    displacements[rank] = total_count;
    total_count += counts[rank];
  }
  std::vector<double> receive(total_count, -1);

  comm->gatherv(send.data(), receive.data(), local_count, counts, displacements);

  if (comm->rank() == 0)
    for (int rank = 0; rank < comm->size(); ++rank)
      for (int element = 0; element < counts[rank]; ++element)
        CHECK(receive[displacements[rank] + element] == 10 * rank + element + 1);
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
      CHECK(values(row, column) == Complex{static_cast<double>(value), static_cast<double>(-value)});
    }
}

TEST_CASE("test_communicate_complex_pointer_reduce_in_place", "[message]")
{
  Communicate* c = OHMMS::Controller;
  std::vector<Complex> values{{(double)(c->rank() + 1), (double)(c->rank() + 2)},
                                           {(double)(c->rank() + 3), (double)(c->rank() + 4)}};

  c->reduce_in_place(values.data(), values.size());

  if (c->rank() == 0)
  {
    double expected_real_0 = 0.0;
    double expected_imag_0 = 0.0;
    double expected_real_1 = 0.0;
    double expected_imag_1 = 0.0;
    for (int i = 0; i < c->size(); i++)
    {
      expected_real_0 += (double)(i + 1);
      expected_imag_0 += (double)(i + 2);
      expected_real_1 += (double)(i + 3);
      expected_imag_1 += (double)(i + 4);
    }
    REQUIRE(values[0] == Complex{expected_real_0, expected_imag_0});
    REQUIRE(values[1] == Complex{expected_real_1, expected_imag_1});
  }
}

TEST_CASE("communicate_collectives_allgather_count", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;

  std::vector<int> send{1 + comm->rank(), 2 + comm->rank(), 3 + comm->rank()};
  std::vector<int> receive(send.size() * comm->size(), 0);
  comm->allgather(send, receive);

  for (int i = 0; i < comm->size(); i++)
  {
    CHECK(receive[i * 3 + 0] == 1 + i);
    CHECK(receive[i * 3 + 1] == 2 + i);
    CHECK(receive[i * 3 + 2] == 3 + i);
  }
}
} // namespace qmcplusplus
