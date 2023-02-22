// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#ifndef BOOST_MPI3_WALL_CLOCK_HPP
#define BOOST_MPI3_WALL_CLOCK_HPP

#include <mpi.h>

#include <mpi3/communicator.hpp>

#include <chrono>

namespace boost {
namespace mpi3 {

inline auto wall_time() { return std::chrono::duration<double>(MPI_Wtime()); }
inline auto wall_tick() { return std::chrono::duration<double>(MPI_Wtick()); }

struct wall_clock {
	using rep        = double;
	using period     = std::ratio<1>;  // one second
	using duration   = std::chrono::duration<double>;
	using time_point = std::chrono::time_point<wall_clock>;

	static time_point now() noexcept { return time_point{wall_time()}; };
	static duration   tick() noexcept { return duration{wall_tick()}; }
};

template<class Duration = std::chrono::duration<double>>
void wall_sleep_for(Duration d) {
	auto then = wall_clock::now();
	// spin now
	while(wall_clock::now() - then < d) {  // NOLINT(altera-unroll-loops) spin loop
	}
}

class wall_timer {
	mpi3::communicator     comm_;
	std::string            title_;
	wall_clock::time_point start_;

 public:
	explicit wall_timer(mpi3::communicator comm, std::string title = "")
	: comm_{std::move(comm)}, title_{std::move(title)}, start_{wall_clock::now()} {}

	wall_timer(wall_timer const&) = delete;
	wall_timer(wall_timer&&)      = delete;

	wall_timer& operator=(wall_timer const&) = delete;
	wall_timer& operator=(wall_timer&&)      = delete;

	~wall_timer() {  // NOLINT(bugprone-exception-escape) TODO(correaa) maybe it should be able to throw
		auto const diff     = wall_clock::now() - start_;
		auto const min      = comm_.min(diff.count());  // cppcheck-suppress unreadVariable ; bug in cppcheck 2.3
		auto const max      = comm_.max(diff.count());  // cppcheck-suppress unreadVariable ; bug in cppcheck 2.3
		auto const total    = (comm_ += diff.count());
		auto const avg      = total / comm_.size();
		auto const speed_up = max / total;
		if(comm_.root()) {
			std::cerr << "# " << title_ << " timing " << min << "[" << avg << "]" << max << " sec, speed up = x" << speed_up << std::endl;
		}
	}
};

}  // end namespace mpi3
}  // end namespace boost

#endif
