// Copyright 2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/window.hpp>
/* This does a transpose-cum-accumulate operation. Uses vector and hvector datatypes (Example 3.32 from MPI 1.1 Standard). Run on 2 processes */

constexpr std::size_t NROWS = 100;
constexpr std::size_t NCOLS = 100;

namespace mpi3 = boost::mpi3;

int main(int argc, char* argv[]) {  // NOLINT(bugprone-exception-escape,readability-function-cognitive-complexity)
	mpi3::environment env(argc, argv);

	auto comm = (env.world() < 2);

	if(comm) {
		assert(comm.size() == 2);

		std::array<std::array<int, NCOLS>, NROWS> A{};

		if(comm.rank() == 0) {
			for(std::size_t i = 0; i < NROWS; ++i) {
				for(std::size_t j = 0; j < NCOLS; ++j) {          // NOLINT(altera-unroll-loops)
					A[i][j] = static_cast<int>((i * NCOLS) + j);  // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
				}
			}

			MPI_Datatype column;  // NOLINT(cppcoreguidelines-init-variables) delayed init
			MPI_Datatype xpose;   // NOLINT(cppcoreguidelines-init-variables) delayed init

			/* create datatype for one column */
			MPI_Type_vector(NROWS, 1, NCOLS, MPI_INT, &column);
			/* create datatype for matrix in column-major order */
			MPI_Type_create_hvector(NCOLS, 1, sizeof(int), column, &xpose);

			MPI_Type_commit(&xpose);

			auto win = comm.create_window<char*>();  // mpi3::communicator

			{
				mpi3::epoch _(win);  // NOLINT(misc-const-correctness)
				MPI_Accumulate(A[0].data(), NROWS * NCOLS, MPI_INT, 1, 0, 1, xpose, MPI_SUM, &win);
			}

			MPI_Type_free(&column);
			MPI_Type_free(&xpose);
		} else {
			assert(comm.rank() == 1);
			for(std::size_t i = 0; i != NROWS; ++i) {
				for(std::size_t j = 0; j != NCOLS; ++j) {         // NOLINT(altera-unroll-loops)
					A[i][j] = static_cast<int>((i * NCOLS) + j);  // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
				}
			}

			auto win = comm.create_window(A[0].data(), NROWS * NCOLS);

			{
				mpi3::epoch const _(win);
			}

			for(std::size_t j = 0; j != NCOLS; ++j) {
				for(std::size_t i = 0; i != NROWS; ++i) {                                    // NOLINT(altera-unroll-loops)
					assert(A[j][i] == static_cast<int>((i * NCOLS) + j + (j * NCOLS) + i));  // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
				}
			}
		}
	}

	return 0;
}
