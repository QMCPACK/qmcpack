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

// =======================================================================
// BCAST API: Scalars and Explicit Specializations
// =======================================================================
TEST_CASE("communicate_collectives_scalar_bcast", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;

  // Scalar double
  double val_d = 0.0;
  if (comm->rank() == 0)
    val_d = 3.0;
  comm->bcast(val_d);
  CHECK(val_d == 3.0);

  // Scalar complex
  Complex val_c{0.0, 0.0};
  if (comm->rank() == 0)
    val_c = {1.2, 3.4};
  comm->bcast(val_c);
  CHECK(val_c == Complex{1.2, 3.4});

  // Explicit bool
  bool val_b = false;
  if (comm->rank() == 0)
    val_b = true;
  comm->bcast(val_b);
  CHECK(val_b == true);

  // Explicit std::string
  std::string val_s = "";
  if (comm->rank() == 0)
    val_s = "qmcpack_mpi_string";
  comm->bcast(val_s);
  CHECK(val_s == "qmcpack_mpi_string");
}

TEST_CASE("communicate_collectives_complex_vector_bcast", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;
  std::vector<Complex> values{{-1, -1}, {-1, -1}};
  if (comm->rank() == 0)
    values = {{1, -2}, {3, -4}};
  comm->bcast(values);
  CHECK(values[0] == Complex{1, -2});
  CHECK(values[1] == Complex{3, -4});
}

TEST_CASE("communicate_collectives_vector_bool_bcast", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;

  std::vector<bool> val_vb(3, false);
  if (comm->rank() == 0)
  {
    val_vb[0] = true;
    val_vb[1] = false;
    val_vb[2] = true;
  }

  comm->bcast(val_vb);

  REQUIRE(val_vb.size() == 3);
  CHECK(val_vb[0] == true);
  CHECK(val_vb[1] == false);
  CHECK(val_vb[2] == true);
}

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

// =======================================================================
// REDUCE API: Scalars
// =======================================================================
TEST_CASE("communicate_collectives_scalar_reduce", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;

  // Scalar int
  int val_i = comm->rank() + 1;
  comm->reduce(val_i);
  if (comm->rank() == 0)
  {
    int expected_i = comm->size() * (comm->size() + 1) / 2;
    CHECK(val_i == expected_i);
  }

  // Scalar double
  double val_d = static_cast<double>(comm->rank() + 1);
  comm->reduce(val_d);
  if (comm->rank() == 0)
  {
    double expected_d = static_cast<double>(comm->size() * (comm->size() + 1) / 2);
    CHECK(val_d == expected_d);
  }

  // Scalar complex
  Complex val_c{static_cast<double>(comm->rank() + 1), static_cast<double>(comm->rank() + 2)};
  comm->reduce(val_c);
  if (comm->rank() == 0)
  {
    double expected_real = static_cast<double>(comm->size() * (comm->size() + 1) / 2);
    double expected_imag = static_cast<double>(comm->size() * (comm->size() + 1) / 2 + comm->size());
    CHECK(val_c == Complex{expected_real, expected_imag});
  }
}

TEST_CASE("communicate_collectives_complex_vector_reduce", "[message][collectives]")
{
  Communicate* c = OHMMS::Controller;
  std::vector<Complex> values{{(double)(c->rank() + 1), (double)(c->rank() + 2)},
                              {(double)(c->rank() + 3), (double)(c->rank() + 4)}};
  c->reduce(values);
  if (c->rank() == 0)
  {
    double expected_real_0 = 0.0, expected_imag_0 = 0.0;
    double expected_real_1 = 0.0, expected_imag_1 = 0.0;
    for (int i = 0; i < c->size(); i++)
    {
      expected_real_0 += (double)(i + 1);
      expected_imag_0 += (double)(i + 2);
      expected_real_1 += (double)(i + 3);
      expected_imag_1 += (double)(i + 4);
    }
    CHECK(values[0] == Complex{expected_real_0, expected_imag_0});
    CHECK(values[1] == Complex{expected_real_1, expected_imag_1});
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

// =======================================================================
// ALLREDUCE API: Scalars
// =======================================================================
TEST_CASE("communicate_collectives_scalar_allreduce", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;

  // Scalar int
  int val_i = comm->rank() + 1;
  comm->allreduce(val_i);
  int expected_i = comm->size() * (comm->size() + 1) / 2;
  CHECK(val_i == expected_i);

  // Scalar double
  double val_d = static_cast<double>(comm->rank() + 1);
  comm->allreduce(val_d);
  double expected_d = static_cast<double>(comm->size() * (comm->size() + 1) / 2);
  CHECK(val_d == expected_d);

  // Scalar complex
  Complex val_c{static_cast<double>(comm->rank() + 1), static_cast<double>(comm->rank() + 2)};
  comm->allreduce(val_c);
  double expected_real = static_cast<double>(comm->size() * (comm->size() + 1) / 2);
  // Sum of (rank + 2):
  // Sum(rank+1) = S. Sum(1) = size. So S + size.
  double expected_imag = static_cast<double>(comm->size() * (comm->size() + 1) / 2 + comm->size());
  CHECK(val_c == Complex{expected_real, expected_imag});
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

TEST_CASE("communicate_collectives_complex_vector_gather", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;
  std::vector<Complex> send{rank_value(comm->rank(), 0), rank_value(comm->rank(), 1)};
  std::vector<Complex> receive;
  if (comm->rank() == 0)
    receive.resize(comm->size() * send.size(), Complex{-1, -1});

  comm->gather(send, receive, 0);

  if (comm->rank() == 0)
  {
    for (int rank = 0; rank < comm->size(); ++rank)
      for (int element = 0; element < static_cast<int>(send.size()); ++element)
        CHECK(receive[rank * send.size() + element] == rank_value(rank, element));
  }
}

TEST_CASE("communicate_collectives_float_vector_gatherv", "[message][collectives]")
{
  Communicate* comm     = OHMMS::Controller;
  const int local_count = comm->rank() + 1;
  std::vector<float> send(local_count);
  for (int element = 0; element < local_count; ++element)
    send[element] = 10.0f * comm->rank() + element + 1.0f;

  std::vector<int> counts(comm->size());
  std::vector<int> displacements(comm->size());
  int total_count = 0;
  for (int rank = 0; rank < comm->size(); ++rank)
  {
    counts[rank]        = rank + 1;
    displacements[rank] = total_count;
    total_count += counts[rank];
  }
  std::vector<float> receive(total_count, -1.0f);

  comm->gatherv(send, receive, counts, displacements);

  if (comm->rank() == 0)
    for (int rank = 0; rank < comm->size(); ++rank)
      for (int element = 0; element < counts[rank]; ++element)
        CHECK(receive[displacements[rank] + element] == 10.0f * rank + element + 1.0f);
}

TEST_CASE("communicate_collectives_double_pointer_gatherv", "[message][collectives]")
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

TEST_CASE("communicate_collectives_float_pointer_gatherv", "[message][collectives]")
{
  Communicate* comm     = OHMMS::Controller;
  const int local_count = comm->rank() + 1;
  std::vector<float> send(local_count);
  for (int element = 0; element < local_count; ++element)
    send[element] = 10.0f * comm->rank() + element + 1.0f;

  std::vector<int> counts(comm->size());
  std::vector<int> displacements(comm->size());
  int total_count = 0;
  for (int rank = 0; rank < comm->size(); ++rank)
  {
    counts[rank]        = rank + 1;
    displacements[rank] = total_count;
    total_count += counts[rank];
  }
  std::vector<float> receive(total_count, -1.0f);

  comm->gatherv(send.data(), receive.data(), local_count, counts, displacements);

  if (comm->rank() == 0)
    for (int rank = 0; rank < comm->size(); ++rank)
      for (int element = 0; element < counts[rank]; ++element)
        CHECK(receive[displacements[rank] + element] == 10.0f * rank + element + 1.0f);
}

TEST_CASE("communicate_collectives_float_pointer_gatherv_in_place", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;
  std::vector<int> counts(comm->size());
  std::vector<int> displacements(comm->size());
  int total_count = 0;
  for (int rank = 0; rank < comm->size(); ++rank)
  {
    counts[rank]        = rank + 1;
    displacements[rank] = total_count;
    total_count += counts[rank];
  }

  std::vector<float> buffer(total_count, -1.0f);

  for (int element = 0; element < counts[comm->rank()]; ++element)
    buffer[displacements[comm->rank()] + element] = 10.0f * comm->rank() + element + 1.0f;

  comm->gatherv_in_place(buffer.data(), qmcplusplus::mpi::get_mpi_datatype(buffer[0]), counts, displacements);

  if (comm->rank() == 0)
    for (int rank = 0; rank < comm->size(); ++rank)
      for (int element = 0; element < counts[rank]; ++element)
        CHECK(buffer[displacements[rank] + element] == 10.0f * rank + element + 1.0f);
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

TEST_CASE("communicate_collectives_complex_vector_allgather", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;
  std::vector<Complex> send{rank_value(comm->rank(), 0), rank_value(comm->rank(), 1)};
  std::vector<Complex> receive(comm->size() * send.size(), Complex{-1, -1});

  comm->allgather(send, receive);

  for (int rank = 0; rank < comm->size(); ++rank)
    for (int element = 0; element < static_cast<int>(send.size()); ++element)
      CHECK(receive[rank * send.size() + element] == rank_value(rank, element));
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

TEST_CASE("communicate_collectives_complex_vector_scatter", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;
  std::vector<Complex> send;
  if (comm->rank() == 0)
  {
    send.resize(comm->size() * 2);
    for (int rank = 0; rank < comm->size(); ++rank)
    {
      send[rank * 2 + 0] = rank_value(rank, 0);
      send[rank * 2 + 1] = rank_value(rank, 1);
    }
  }
  std::vector<Complex> receive(2, Complex{-1, -1});

  comm->scatter(send, receive, 0);

  CHECK(receive[0] == rank_value(comm->rank(), 0));
  CHECK(receive[1] == rank_value(comm->rank(), 1));
}

TEST_CASE("communicate_collectives_float_vector_scatterv", "[message][collectives]")
{
  Communicate* comm = OHMMS::Controller;
  std::vector<int> counts(comm->size());
  std::vector<int> displacements(comm->size());
  int total_count = 0;
  for (int rank = 0; rank < comm->size(); ++rank)
  {
    counts[rank]        = rank + 1;
    displacements[rank] = total_count;
    total_count += counts[rank];
  }

  std::vector<float> send;
  if (comm->rank() == 0)
  {
    send.resize(total_count, -1.0f);
    for (int rank = 0; rank < comm->size(); ++rank)
      for (int element = 0; element < counts[rank]; ++element)
        send[displacements[rank] + element] = 10.0f * rank + element + 1.0f;
  }

  const int local_count = comm->rank() + 1;
  std::vector<float> receive(local_count, -1.0f);

  comm->scatterv(send, receive, counts, displacements, 0);

  for (int element = 0; element < local_count; ++element)
    CHECK(receive[element] == 10.0f * comm->rank() + element + 1.0f);
}

} // namespace qmcplusplus
