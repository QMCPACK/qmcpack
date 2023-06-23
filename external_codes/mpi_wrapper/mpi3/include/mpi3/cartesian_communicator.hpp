//  -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#ifndef BOOST_MPI3_CARTESIAN_COMMUNICATOR_HPP
#define BOOST_MPI3_CARTESIAN_COMMUNICATOR_HPP

#include <mpi3/communicator.hpp>
#include <mpi3/process.hpp>

#include <mpi3/detail/call.hpp>

namespace boost::mpi3 {

using dimensionality_type = int;

static constexpr dimensionality_type dynamic_extent = -1;

template<dimensionality_type D = dynamic_extent> struct cartesian_communicator;

template<>
struct cartesian_communicator<dynamic_extent> : communicator {

	cartesian_communicator() = default;

	cartesian_communicator(cartesian_communicator const&) = delete;
	cartesian_communicator(cartesian_communicator&&)      = default;
	// vvv---  this is an unusual "duplicate" constructor
	cartesian_communicator(cartesian_communicator& other) : communicator{other} {}  // NOLINT(hicpp-use-equals-default,modernize-use-equals-default) cannot be defaulted because bug in nvcc 11

	template<class Shape, class Period>
	cartesian_communicator(communicator& comm_old, Shape const& s, Period const& p) {
		assert(s.size() == p.size());
		using dimensionality_type = int;
		MPI_(Cart_create)(comm_old.get(), static_cast<dimensionality_type>(s.size()), s.data(), p.data(), /*reorder*/ true, &impl_);
		//	assert(impl_ != MPI_COMM_NULL); // null communicator is a valid outcome
		// TODO(correaa) try with mpich, WAS: there is an bug in mpich, in which if the remaining dim are none then the communicator is not well defined.
	}

	template<class Shape>
	cartesian_communicator(communicator& old, Shape const& s)
	: cartesian_communicator{old, s, std::vector<int>(s.size(), true)} {}

	cartesian_communicator(communicator& comm_old, std::initializer_list<int> shape)
	: cartesian_communicator(comm_old, std::vector<int>(shape)) {}

	cartesian_communicator(communicator& comm_old, std::initializer_list<int> shape, std::initializer_list<int> period)
	: cartesian_communicator(comm_old, std::vector<int>(shape), std::vector<int>(period)) {}

	[[deprecated("use dimensionality() instead of dimension")]] int dimension() const {
		int ret;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		MPI_Cartdim_get(impl_, &ret);
		return ret;
	}

	cartesian_communicator& operator=(cartesian_communicator const&) = delete;
	cartesian_communicator& operator=(cartesian_communicator&&)      = default;
	// vvv nvcc 11 workaround, needs explicit definition of duplicate assigment
	[[deprecated]] cartesian_communicator& operator=(cartesian_communicator& other) {  // NOLINT(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator) "duplicate" assignment
		if(this == std::addressof(other)) {
			return *this;
		}  // lints cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator
		communicator::operator=(other);
		return *this;
	}

	~cartesian_communicator() = default;

	int dimensionality() const {
		int ret;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		MPI_(Cartdim_get)(impl_, &ret);
		return ret;
	}

	std::vector<int> coordinates() const {
		std::vector<int> ret(static_cast<std::vector<int>::size_type>(dimensionality()));
		MPI_(Cart_coords)(impl_, rank(), dimensionality(), ret.data());
		return ret;
	}

	auto topology() const {
		auto const maxdims = static_cast<std::size_t>(dimensionality());  // TODO(correaa) use safe cast
		class topology_t {
			std::vector<int> dimensions_;
			std::vector<int> periods_;
			std::vector<int> coordinates_;
			friend mpi3::cartesian_communicator<dynamic_extent>;

		 public:
			explicit topology_t(std::size_t n) : dimensions_(n), periods_(n), coordinates_(n) {}

			auto const& dimensions() const { return dimensions_; }
			auto const& periods() const { return periods_; }
			auto const& coordinates() const { return coordinates_; }
		} ret(maxdims);

		MPI_(Cart_get)(impl_, static_cast<dimensionality_type>(maxdims), ret.dimensions_.data(), ret.periods_.data(), ret.coordinates_.data());

		assert(ret.coordinates() == coordinates());
		return ret;
	}

	std::vector<int>  shape() const { return topology().dimensions(); }
	std::vector<bool> periods() const {
		auto ps = topology().periods();
		return {ps.begin(), ps.end()};
	}
	auto num_elements() const { return size(); }

	template<class Coord>
	auto operator()(Coord const& coord) {
		int rank = -1;
		MPI_(Cart_rank)(impl_, coord.data(), &rank);
		return (*this)[rank];
		//	return operator[](rank);
	}
	// int MPI_Cart_map not implemented
	cartesian_communicator sub_aux(std::vector<int> const& remain_dims) {
		assert(static_cast<dimensionality_type>(remain_dims.size()) == dimensionality());
		cartesian_communicator ret;
		MPI_(Cart_sub)(impl_, remain_dims.data(), &ret.impl_);
		return ret;
	}

	template<class RemainDim = std::initializer_list<bool>>
	cartesian_communicator sub(RemainDim const& remain_dims) {
		return sub_aux(std::vector<int>(remain_dims.begin(), remain_dims.end()));
	}
	cartesian_communicator sub() {
		assert(dimensionality() > 1);
		std::vector<int> remain(static_cast<std::size_t>(dimensionality()), 1 /*true*/);
		remain[0] = 0 /*false*/;
		return sub_aux(remain);
	}
};

enum fill_t {
	fill = 0,
	_    = 0
};

struct circular_communicator;

template<dimensionality_type D>
struct cartesian_communicator : cartesian_communicator<> {
	cartesian_communicator() = default;

	cartesian_communicator(cartesian_communicator& other) : cartesian_communicator<>{other} {}
	cartesian_communicator(cartesian_communicator const&)     = delete;
	cartesian_communicator(cartesian_communicator&&) noexcept = default;

	~cartesian_communicator() = default;

	static std::array<int, D> division(int nnodes, std::array<int, D> suggest = {}) {
		MPI_(Dims_create)(nnodes, D, suggest.data());
		return suggest;
	}
	constexpr static dimensionality_type dimensionality = D;

	explicit cartesian_communicator(
		communicator&       other,
		std::array<int, D>  dims    = {},
		std::array<bool, D> periods = std::apply([](auto... e) { return std::array{(static_cast<void>(e), true)...}; }, std::array<int, D>{})
	) try : cartesian_communicator
		<>{other, division(other.size(), dims), std::apply([](auto... e) { return std::array<int, D>{e...}; }, periods)} {}
	catch(std::runtime_error& e) {
		std::ostringstream ss;
		std::copy(dims.begin(), dims.end(), std::ostream_iterator<int>{ss, " "});
		throw std::runtime_error{"cannot create cartesian communicator with constrains " + ss.str() + " from communicator of size " + std::to_string(other.size()) + " because " + e.what()};
	}

	auto topology() const {
		struct topology_t {
			std::array<int, dimensionality> dimensions, periods, coordinates;
		} ret = {};
		MPI_(Cart_get)(
			impl_, dimensionality,
			ret.dimensions.data(), ret.periods.data(), ret.coordinates.data()
		);
		return ret;
	}

	constexpr auto dimensions() const { return topology().dimensions; }

	cartesian_communicator& operator=(cartesian_communicator const&)     = delete;
	cartesian_communicator& operator=(cartesian_communicator&&) noexcept = default;
	// vvv  nvcc 11 workaround, needs explicit definition of duplicate assigment
	[[deprecated]] cartesian_communicator& operator=(cartesian_communicator& other) {  // NOLINT(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator) duplicate assignment
		if(this == std::addressof(other)) {
			return *this;
		}  // lints cert-oop54-cpp
		cartesian_communicator<>::operator=(other);  // NOLINT(clang-diagnostic-deprecated-declarations)
		return *this;
	}

	cartesian_communicator<1> axis(int d) {  // TODO(correaa) return a circular_communicator
		assert(d >= 0);
		assert(d < D);
		cartesian_communicator<1> ret;
		std::array<int, D>        remains{};
		remains.fill(false);
		remains[static_cast<std::size_t>(d)] = true;  // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
		MPI_(Cart_sub)(impl_, remains.data(), &ret.get());
		return ret;
	}

	template<int DD> auto axis() -> circular_communicator;

	auto plane(int d1, int d2) {
		assert(d1 >= 0 and d2 >= 0);

		auto const d1s = static_cast<std::size_t>(d1);
		auto const d2s = static_cast<std::size_t>(d2);

		assert(d2s < D and d1s < d2s);

		std::array<int, D> remains{};
		remains.at(d1s) = true;
		remains.at(d2s) = true;

		cartesian_communicator<2> ret;
		MPI_(Cart_sub)(impl_, remains.data(), &ret.get());
		return ret;
	}

	template<int D1 = 0, int D2 = 1>
	auto plane() {
		static_assert(D1 < D2);
		std::array<int, D> remains{};
		remains.fill(false);
		std::get<D1>(remains) = true;
		std::get<D2>(remains) = true;

		cartesian_communicator<2> ret;
		MPI_(Cart_sub)(impl_, remains.data(), &ret.get());
		return ret;
	}

	template<int Direction> auto shift(int displacement) const {
		std::pair<int, int> source_dest;
		MPI_(Cart_shift)(impl_, Direction, displacement, &source_dest.first, &source_dest.second);
		return source_dest;
	}

	template<int Direction, typename... As>
	auto send_receive_shift(As... as, int displacement = 1) {
		send_receive(as..., shift<Direction>(displacement));
	}

	using coordinates_type = std::array<int, D>;

	using cartesian_communicator<>::rank;
	auto rank(coordinates_type cs) const -> int {
		auto const ps = periods();
		auto const s  = shape();
		for(std::size_t i = 0; i != D; ++i) {  // NOLINT(altera-unroll-loops) TODO(correaa) use algorithm
			if(ps[i] == false) {
				assert(cs[i] >= 0);
				assert(cs[i] < s[i]);
			}
		}
		return MPI_(Cart_rank)(impl_, cs.data());
	}
	auto coordinates(int r) const -> coordinates_type {
		coordinates_type ret;
		MPI_(Cart_coords)(impl_, r, D, ret.data());
		return ret;
	}
	auto coordinates() const -> coordinates_type { return coordinates(rank()); }
	template<class... Indices>
	auto operator()(Indices... idx) { return (*this)[rank(std::array<int, D>{idx...})]; }

	cartesian_communicator<D - 1> hyperplane(int d) {
		static_assert(D != 0, "hyperplane not possible for 0D communicators");
#if defined(MPICH_VERSION)
		static_assert(D != 1, "hyperplane not possible for 1D communicators");  // they work in openMPI but do not work in MPICH
#endif
		assert(d >= 0);
		assert(d < D);

		cartesian_communicator<D - 1> ret;
		std::array<int, D>            remains{};
		remains.fill(true);

		remains[static_cast<std::size_t>(d)] = false;  // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
		MPI_(Cart_sub)(impl_, remains.data(), &ret.get());
		return ret;
	}
};

struct circular_communicator : cartesian_communicator<1> {
	circular_communicator() = default;

	circular_communicator(circular_communicator& other) : cartesian_communicator<1>{other} {}  // NOLINT(hicpp-use-equals-default,modernize-use-equals-default) icpc needs this definition explicitly
	circular_communicator(circular_communicator const&)     = delete;
	circular_communicator(circular_communicator&&) noexcept = default;

	~circular_communicator() = default;

	explicit circular_communicator(communicator& other)
	: cartesian_communicator<1>{other, {}, {true}} {}

	auto operator=(cartesian_communicator const&) -> circular_communicator& = delete;
	auto operator=(circular_communicator&& other) noexcept -> circular_communicator& {
		cartesian_communicator<1>::operator=(std::move(other));
		return *this;
	}
	// vvv  nvcc 11 workaround, needs explicit definition of duplicate assigment
	[[deprecated]] auto operator=(circular_communicator& other) -> circular_communicator& {  // NOLINT(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator) duplicate assignment
		if(this == std::addressof(other)) {
			return *this;
		}  // lints cert-oop54-cpp
		cartesian_communicator<1>::operator=(other);  // NOLINT(clang-diagnostic-deprecated-declarations)
		return *this;
	}

	auto coordinate() const { return std::get<0>(this->coordinates()); }
	auto coordinate(int rank) const { return std::get<0>(this->coordinates(rank)); }

	using cartesian_communicator<1>::rank;
	auto rank(int coordinate) const { return cartesian_communicator<1>::rank({coordinate}); }

	template<typename... As>
	auto rotate(As... as, int displacement) { return this->send_receive(as..., this->shift<0>(-displacement)); }
	template<typename... As>
	auto rotate(As... as) { return this->send_receive(as..., this->shift<0>(-1)); }

	template<typename... As>
	auto unrotate(As... as, int displacement) { return this->send_receive(as..., this->shift<0>(+displacement)); }
	template<typename... As>
	auto unrotate(As... as) { return this->send_receive(as..., this->shift<0>(+1)); }
};

template<int Dimension> using torus = cartesian_communicator<Dimension>;

using ring = circular_communicator;

template<int D> template<int DD>
auto cartesian_communicator<D>::axis() -> circular_communicator {
	circular_communicator ret;
	std::array<int, D>    remains{};
	remains.fill(false);
	std::get<DD>(remains) = true;
	MPI_(Cart_sub)(impl_, remains.data(), &ret.get());
	return ret;
}

}  // end namespace boost::mpi3
#endif
