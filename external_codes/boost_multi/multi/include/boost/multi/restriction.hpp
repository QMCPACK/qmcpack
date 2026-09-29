// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_RESTRICTION_HPP
#define BOOST_MULTI_RESTRICTION_HPP

#include "boost/multi/detail/layout.hpp"  // IWYU pragma: export
#include "boost/multi/utility.hpp"

#ifdef __NVCC__
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
#endif
// NOLINTNEXTLINE(bugprone-macro-parentheses)
#define BOOST_MULTI_HD_LAMBDA(BodY)     \
	{                                   \
		return [=] BOOST_MULTI_HD BodY; \
	}                                   \
	()

namespace boost::multi::detail {

#ifdef __NVCC__
template<class Fun>
struct device : Fun {
	using Fun::operator();
};
#endif

template<class Fun> struct function_system {
	using type = void;
};

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

template<class T>
constexpr auto make_restriction(std::initializer_list<T> const& ilist) {
	return [ilist](multi::index idx0) -> decltype(auto) { return ilist.begin()[idx0]; } ^ multi::extents(ilist);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

template<class T>
constexpr auto make_restriction(std::initializer_list<std::initializer_list<T>> const& ilist) {
	return [ilist](multi::index idx0, multi::index idx1) -> decltype(auto) { return ilist.begin()[idx0].begin()[idx1]; } ^ multi::extents(ilist);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

template<class T>
constexpr auto make_restriction(std::initializer_list<std::initializer_list<std::initializer_list<T>>> const& ilist) {
	return [ilist](multi::index idx0, multi::index idx1, multi::index idx2) -> decltype(auto) { return ilist.begin()[idx0].begin()[idx1].begin()[idx2]; } ^ multi::extents(ilist);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

template<class T, dimensionality_type D>
struct init_list {
	using type = std::initializer_list<typename init_list<T, D - 1>::type>;
};

template<class T>
struct init_list<T, 1> {
	using type = std::initializer_list<T>;
};

template<class T>
struct init_list<T, 0> {
	using type = T;
};

template<class T, dimensionality_type D>
using init_list_t = typename init_list<T, D>::type;

}  // namespace boost::multi::detail

namespace boost::multi {

namespace detail {
template<class T, dimensionality_type D>
using restriction_idl = decltype(detail::make_restriction(std::declval<detail::init_list_t<T, D>>()));

template<class T, dimensionality_type D>
class initializer_array : public restriction_idl<T, D> {
	using base_ = restriction_idl<T, D>;
	// detail::init_list_t<T, D> ild_;
 public:
	// cppcheck-suppress noExplicitConstructor ;
	constexpr initializer_array(detail::init_list_t<T, D> ild)  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,cppcoreguidelines-explicit-constructor,misc-explicit-constructor)
	: base_(detail::make_restriction(ild)) {}
};

template<class Proj>
struct bind_transposed_t {
	Proj proj_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

	template<class T1, class T2, class... Ts>
	BOOST_MULTI_HD constexpr auto operator()(T1 row, T2 col, Ts... rest) const noexcept -> decltype(auto) /*element*/ { return proj_(col, row, rest...); }
};

template<class Proj>
struct bind_front_t {
	multi::index idx_;
	Proj         proj_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

	// bind_front_t(multi::index idx, Proj& proj) : idx_{idx}, proj_{proj} {}
	template<class... Args>
	BOOST_MULTI_HD constexpr auto operator()(Args&&... rest) const& noexcept -> decltype(auto) { return proj_(idx_, std::forward<Args>(rest)...); }

	template<class... Args>
	BOOST_MULTI_HD constexpr auto operator()(Args&&... rest) & noexcept -> decltype(auto) { return proj_(idx_, std::forward<Args>(rest)...); }

	template<class... Args>
	BOOST_MULTI_HD constexpr auto operator()(Args&&... rest) && noexcept -> decltype(auto) { return proj_(idx_, std::forward<Args>(rest)...); }
};
}  // end namespace detail

template<dimensionality_type D, class Proj>
class restriction;

namespace detail {
template<class F, class TupleLike>
struct invoke_result_from_tuple;

template<class F, template<class...> class TupleLike, class... Args>
struct invoke_result_from_tuple<F, TupleLike<Args...>> {
	using type =
		std::invoke_result_t<F, std::decay_t<Args>...>;
};
}  // end namespace detail

namespace detail {
/// Iterator for a restriction array in the leading dimension
template<dimensionality_type D, class Proj>
class restriction_iterator {
	typename extents_t<D>::iterator it_;
	Proj const*                     Pproj_;

	restriction_iterator(typename extents_t<D>::iterator iter, Proj const* Pproj) : it_{iter}, Pproj_{Pproj} {}

	template<dimensionality_type, class>
	friend class restriction;

	template<class Fun, class... Args>
	static BOOST_MULTI_HD constexpr auto apply_(Fun&& fun, Args&&... args) {  // NOLINT(readability-identifier-naming) TODO(correaa) move this to ::detail
		using std::apply;
		return apply(std::forward<Fun>(fun), std::forward<Args>(args)...);
	}

 public:
	constexpr restriction_iterator() = default;  // cppcheck-suppress uninitMemberVar ; partially formed

	restriction_iterator(restriction_iterator const& other) = default;
	restriction_iterator(restriction_iterator&&) noexcept   = default;

	auto operator=(restriction_iterator&&) noexcept -> restriction_iterator& = default;
	auto operator=(restriction_iterator const&) -> restriction_iterator&     = default;

	~restriction_iterator() = default;

	/// Tag the execution system
	using system = typename multi::detail::function_system<Proj>::type;

	/// Signed integer type for iterator arithmetic
	using difference_type   = std::ptrdiff_t;
	/// Value type resulting from binding the leading index, or, for `D == 1` the element type
	using value_type        = std::conditional_t<(D != 1), restriction<D - 1, detail::bind_front_t<Proj>>, decltype(apply_(std::declval<Proj>(), std::declval<typename extents_t<D>::element>()))>;
	/// Pointer type (`void` because there is no memory behind)
	using pointer           = void;
	/// Reference type resulting from binding the leading index, or, for `D == 1` the element type
	using reference         = std::conditional_t<(D != 1), restriction<D - 1, detail::bind_front_t<Proj>>, decltype(apply_(std::declval<Proj&>(), std::declval<typename extents_t<D>::element>()))>;
	/// Iterator category is random-access
	using iterator_category = std::random_access_iterator_tag;

	constexpr auto operator++() -> auto& {
		++it_;
		return *this;
	}
	constexpr auto operator--() -> auto& {
		--it_;
		return *this;
	}

	/// Addition assignment operators
	constexpr auto operator+=(difference_type diff) -> auto& {
		it_ += diff;
		return *this;
	}

	/// Substract assignment operators
	constexpr auto operator-=(difference_type diff) -> auto& {
		it_ -= diff;
		return *this;
	}

	constexpr auto operator++(int) {
		restriction_iterator ret{*this};
		++(*this);
		return ret;
	}
	constexpr auto operator--(int) {
		restriction_iterator ret{*this};
		--(*this);
		return ret;
	}

	friend constexpr auto operator-(restriction_iterator const& self, restriction_iterator const& other) { return self.it_ - other.it_; }
	friend constexpr auto operator+(restriction_iterator const& self, difference_type n) {
		restriction_iterator ret{self};
		return ret += n;
	}
	friend constexpr auto operator-(restriction_iterator const& self, difference_type n) {
		restriction_iterator ret{self};
		return ret -= n;
	}

	friend constexpr auto operator+(difference_type n, restriction_iterator const& self) { return self + n; }

	friend constexpr auto operator==(restriction_iterator const& self, restriction_iterator const& other) noexcept -> bool { return self.it_ == other.it_; }
	friend constexpr auto operator!=(restriction_iterator const& self, restriction_iterator const& other) noexcept -> bool { return self.it_ != other.it_; }

	friend auto operator<=(restriction_iterator const& self, restriction_iterator const& other) noexcept -> bool { return self.it_ <= other.it_; }
	friend auto operator<(restriction_iterator const& self, restriction_iterator const& other) noexcept -> bool { return self.it_ < other.it_; }
	friend auto operator>(restriction_iterator const& self, restriction_iterator const& other) noexcept -> bool { return self.it_ > other.it_; }
	friend auto operator>=(restriction_iterator const& self, restriction_iterator const& other) noexcept -> bool { return self.it_ > other.it_; }

	/// Dereference operator
	BOOST_MULTI_HD constexpr auto operator*() const -> decltype(auto) {
		if constexpr(D != 1) {
			using std::get;
			return restriction<D - 1, detail::bind_front_t<Proj>>(extents_t<D - 1>(it_->tail()), detail::bind_front_t<Proj>(get<0>(*it_), *Pproj_));
		} else {
			using std::get;
			return (*Pproj_)(get<0>(*it_));
		}
	}

	/// Subscript operator
	BOOST_MULTI_HD auto operator[](difference_type diff) const { return *((*this) + diff); }  // TODO(correaa) use ra_iterator_facade
};
}  // end namespace detail

namespace detail {
/// A random-access iterator to go through all the elements of a restriction
template<dimensionality_type D, class Proj>
class restriction_elements_iterator : ra_iterable<restriction_elements_iterator<D, Proj>> {
	typename extents_t<D>::elements_t::iterator it_;
	BOOST_MULTI_NO_UNIQUE_ADDRESS Proj          proj_;

 public:
	restriction_elements_iterator() = default;
	restriction_elements_iterator(typename extents_t<D>::elements_t::iterator iter, Proj proj) : it_{iter}, proj_{std::move(proj)} {}

	/// Increment operator
	auto operator++() -> auto& {
		++this->it_;
		return *this;
	}

	/// Decrement operator
	auto operator--() -> auto& {
		--this->it_;
		return *this;
	}

	/// Addition assignment operator
	constexpr auto operator+=(difference_type diff) -> auto& {
		this->it_ += diff;
		return *this;
	}

	/// Subtraction assignment operator
	constexpr auto operator-=(difference_type diff) -> auto& {
		this->it_ -= diff;
		return *this;
	}

	/// Difference operator
	friend constexpr auto operator-(restriction_elements_iterator const& self, restriction_elements_iterator const& other) { return self.it_ - other.it_; }

	/// Tag to deduce system of the function (compute space), for compatibility with CUDA
	using system = typename detail::function_system<std::decay_t<Proj>>::type;

	/// Signed integer type for iterator arithmetic
	using difference_type   = std::ptrdiff_t;  // elements_t::difference_type;
	/// Value (element) type, e.g. `T`
	using value_type        = difference_type;
	/// Pointer type, `void` because there is no memory behind
	using pointer           = void;
	/// Referece type, also `T` (not `T&` because there is no memory behind)
	using reference         = value_type;
	/// Iterator category is Random-Access
	using iterator_category = std::random_access_iterator_tag;

	friend auto operator==(restriction_elements_iterator const& self, restriction_elements_iterator const& other) -> bool { return self.it_ == other.it_; }
	friend auto operator!=(restriction_elements_iterator const& self, restriction_elements_iterator const& other) -> bool { return self.it_ != other.it_; }

	friend auto operator<=(restriction_elements_iterator const& self, restriction_elements_iterator const& other) -> bool { return self.it_ <= other.it_; }
	friend auto operator<(restriction_elements_iterator const& self, restriction_elements_iterator const& other) -> bool { return self.it_ < other.it_; }

	/// Dereference operator
	BOOST_MULTI_HD constexpr auto operator*() const -> decltype(auto) {
		using std::apply;
		return apply(proj_, *this->it_);
	}

	/// Subscript operator, same as `*(*this + diff)`
	BOOST_MULTI_HD constexpr auto operator[](difference_type diff) const -> decltype(auto) { return *((*this) + diff); }  // TODO(correaa) use ra_iterator_facade
};
}  // end namespace detail

namespace detail {
/// A random-access range for all the elements of a restriction
template<dimensionality_type D, class Proj>
class restriction_elements_t {
	typename extents_t<D>::elements_t elems_;
	Proj                              proj_;

 public:
	restriction_elements_t(typename extents_t<D>::elements_t elems, Proj proj) : elems_{elems}, proj_{std::move(proj)} {}

	/// Subscript operator
	BOOST_MULTI_HD constexpr auto operator[](index idx) const -> decltype(auto) {
		using std::apply;
		return apply(proj_, elems_[idx]);
	}

	/// Signed integer type to do iterator arithmetic
	using difference_type = std::ptrdiff_t;

	/// A random-access iterator to access the elements of the restriction
	using iterator  = detail::restriction_elements_iterator<D, Proj>;
	/// Integer type to store the size of the range
	using size_type = typename extents_t<D>::size_type;

	/// returns an output random-iterator to the beggining of the elements range
	auto begin() const -> iterator { return {elems_.begin(), proj_}; }
	/// returns an output random-iterator to the end of the elements range
	auto end() const -> iterator { return {elems_.end(), proj_}; }

	/// returns the size of the elements range
	auto size() const noexcept -> size_type { return elems_.size(); }
};
}  // end namespace detail

/// An array interface for a function of `D` integer arguments restricted to an certain Cartesian grid (extents), elements are generated lazily
template<dimensionality_type D, class Proj>
class restriction : std::conditional_t<std::is_reference_v<Proj>, detail::non_copyable_base, detail::copyable_base> {
	extents_t<D> xs_;
	Proj         proj_;

	using system = typename multi::detail::function_system<std::decay_t<Proj>>::type;

	template<class Fun, class Tup>
	static BOOST_MULTI_HD constexpr auto std_apply_(Fun&& fun, Tup&& tup) -> decltype(auto) {
		using std::apply;
		return apply(std::forward<Fun>(fun), std::forward<Tup>(tup));
	}

 public:
	/// Integer constant describing the dimensionality of the restriction array
	static constexpr dimensionality_type dimensionality = D;
	/// Integer constant describing the dimensionality of the restriction array
	static constexpr dimensionality_type rank_v         = D;

	/// Signed integer type for index arithmetic
	using difference_type = typename extents_t<D>::difference_type;
	/// Integer type to store an index in the leading dimension
	using index           = typename extents_t<D>::index;
	/// Integer type to store the size in the leading dimension
	using size_type       = typename extents_t<D>::size_type;

	/// Tuple to store a set of `D` integer indices to obtain an element
	using indices_type [[deprecated("use indices")]] = typename extents_t<D>::indices_type;

	/// Tuple to store a set of `D` integer indices to obtain an element
	using indices = typename extents_t<D>::element;

	/// constructs a `D`-dimensional restriction array from an extensions argument and a function of `D` arguments
	BOOST_MULTI_HD constexpr restriction(extents_t<D> const& exts, Proj proj) : xs_{exts}, proj_{std::move(proj)} {}

 private:
	using element_ref_  = typename detail::invoke_result_from_tuple<Proj, typename extents_t<D>::element>::type;
	using element_type_ = std::decay_t<typename detail::invoke_result_from_tuple<Proj, typename extents_t<D>::element>::type>;

 public:
	/// Element type of the restriction array
	using element                                    = std::decay_t<typename detail::invoke_result_from_tuple<Proj, typename extents_t<D>::element>::type>;
	/// Element type of the restriction array
	using element_type [[deprecated("use element")]] = element;

	/// Type generated by binding the first index in the array (leading dimension), for `D == 1` this is an element type
	using value_type = std::conditional_t<
		(D == 1),
		element,
		array<element, D - 1>>;

	/// Type generated by binding the first index in the array (leading dimension), for `D == 1` this is an element type
	using reference = std::conditional_t<
		(D != 1),                                               //
		restriction<D - 1, detail::bind_front_t<Proj const&>>,  //
		element                                                 // element_ref_
		>;

#if defined(__cpp_multidimensional_subscript) && (__cpp_multidimensional_subscript >= 202110L)
	template<class... Indices>
	BOOST_MULTI_HD constexpr auto operator[](index idx, Indices... rest) const {
		return operator[](idx)[rest...];
	}
	BOOST_MULTI_HD constexpr auto operator[]() const -> decltype(auto) { return proj_(); }
#endif

	BOOST_MULTI_HD constexpr auto operator[](index idx) && -> decltype(auto) {
		// assert( extent().contains(idx) );
		if constexpr(D != 1) {
			return detail::bind_front_t<Proj>{idx, std::move(proj_)} ^ extents_t<D - 1>(xs_.base().tail());
		} else {
			return proj_(idx);
		}
	}

	BOOST_MULTI_HD constexpr auto operator[](index idx) & -> decltype(auto) {
		// assert( extent().contains(idx) );
		if constexpr(D != 1) {
			return detail::bind_front_t<Proj>{idx, proj_} ^ extents_t<D - 1>(xs_.base().tail());
		} else {
			return proj_(idx);
		}
	}

	BOOST_MULTI_HD constexpr auto operator[](index idx) const& -> decltype(auto) {
		// assert( extent().contains(idx) );
		if constexpr(D != 1) {
			return detail::bind_front_t<Proj const&>{idx, proj_} ^ extents_t<D - 1>(xs_.base().tail());
		} else {
			return proj_(idx);
		}
	}

	/// returns a materialized owning array of the restriction (allocates, O(N) operation)
	constexpr auto operator+() const { return multi::array<element, D>(*this); }

	/// yields a transposed view of the restriction array (first two indices exchanged). (Only for `D >= 2`)
	BOOST_MULTI_HD constexpr auto transposed() && {
		return detail::bind_transposed_t<Proj>{std::move(proj_)} ^ layout_t<D>(extents()).transpose().extents();
		// return [proj = proj_](auto i, auto j, auto... rest) { return proj(j, i, rest...); } ^ layout_t<D>(extents()).transpose().extents();
	}

	BOOST_MULTI_HD constexpr auto transposed() const& -> restriction<D, detail::bind_transposed_t<Proj const&>> {
		return detail::bind_transposed_t<Proj const&>{proj_} ^ layout_t<D>(extents()).transpose().extents();
		// return [proj = proj_](auto i, auto j, auto... rest) { return proj(j, i, rest...); } ^ layout_t<D>(extents()).transpose().extents();
	}

 private:
	struct bind_diagonal_t {
		Proj proj_;
		template<class T1, class... Ts>
		BOOST_MULTI_HD constexpr auto operator()(T1 row, Ts... rest) const noexcept -> element { return proj_(row, row, rest...); }
	};

 public:
	/// yields an array restriction of lower dimensionality that holds the diagonal (in the first two dimensions), array must be square (first two sizes equal). (Only for `D >= 2`).
	BOOST_MULTI_HD constexpr auto diagonal() const -> restriction<D - 1, bind_diagonal_t> {
		static_assert(D > 1);
		using std::get;  // needed for C++17
		return bind_diagonal_t{proj_} ^ std::min(get<0>(sizes()), get<1>(sizes())) * extents().sub().sub();
		// return [proj = proj_](auto i, auto j, auto... rest) { return proj(j, i, rest...); } ^ layout_t<D>(extents()).transpose().extents();
	}

	BOOST_MULTI_HD constexpr auto operator~() && { return std::move(*this).transposed(); }
	BOOST_MULTI_HD constexpr auto operator~() const& { return transposed(); }

 private:
	struct bind_repeat_t {
		Proj proj_;
		template<class... Ts>
		BOOST_MULTI_HD constexpr auto operator()(multi::index /*unused*/, Ts... rest) const noexcept -> element { return proj_(rest...); }
	};

 public:
	/// yields a new restriction of higher dimension by repeating `count` times the restriction as a new leading dimension
	BOOST_MULTI_HD auto repeated(size_type count) const -> restriction<D + 1, bind_repeat_t> {
		return bind_repeat_t{proj_} ^ count * extents();
	}

 private:
	struct bind_partitioned_t {
		Proj      proj_;
		size_type block_size_;
		template<class T1, class T2, class... Ts>
		BOOST_MULTI_HD constexpr auto operator()(T1 row, T2 col, Ts... rest) const noexcept -> element { return proj_(row * block_size_ + col, rest...); }
	};

 public:
	/// yields a restriction of higher dimension by splitting the leading dimension into equal‐sized partitions of size `count` (`count` must divide `size()`)
	BOOST_MULTI_HD constexpr auto partitioned(size_type count) const noexcept -> restriction<D + 1, bind_partitioned_t> {
		return bind_partitioned_t{proj_, size() / count} ^ layout_t<D>(extents()).partition(count).extents();
	}

 private:
	struct bind_reversed_t {
		Proj      proj_;
		size_type size_m1;
		template<class T1, class... Ts>
		BOOST_MULTI_HD constexpr auto operator()(T1 row, Ts... rest) const noexcept -> element { return proj_(size_m1 - row, rest...); }
	};

	BOOST_MULTI_HD constexpr auto reversed() const { return bind_reversed_t{proj_, size() - 1} ^ extents(); }  // NOLINT(readability-identifier-naming)

	struct bind_rotated_t {
		Proj      proj_;
		size_type size_;
		template<class T1, class... Ts>
		BOOST_MULTI_HD constexpr auto operator()(T1 row, Ts... rest) const noexcept { return proj_(rest..., row); }
	};

	struct bind_unrotated_t {
		Proj      proj_;
		size_type size_;
		template<class... Ts, class TN>
		BOOST_MULTI_HD constexpr auto operator()(Ts... rest, TN col) const noexcept { return proj_(col, rest...); }
	};

 public:
	/// yields a restruction that has its indices rotated (to the left)
	BOOST_MULTI_HD constexpr auto rotated() const { return bind_rotated_t{proj_, size()} ^ extents().rotate(); }
	/// yields a restruction that has its indices unrotated (rotated to the right)
	BOOST_MULTI_HD constexpr auto unrotated() const { return bind_unrotated_t{proj_, size()} ^ extents().unrotate(); }

 private:
	template<class Proj2>
	struct bind_element_transformed_t {
		Proj  proj_;
		Proj2 proj2_;
		template<class... Ts>
		BOOST_MULTI_HD constexpr auto operator()(Ts... rest) const noexcept -> element { return proj2_(proj_(rest...)); }
	};

 public:
	/// An array interface with of the same dimensionality to which each element is applied a projection function (unary function) (i.e. each element is obtained by the composition of the restricted function and the projection)
	template<class ElementProj>
	BOOST_MULTI_HD auto element_transformed(ElementProj element_proj) const -> restriction<D, bind_element_transformed_t<ElementProj>> {
		return bind_element_transformed_t<ElementProj>{proj_, element_proj} ^ extents();
	}

 private:
	template<class Proj2>
	class bind_transform_t {
		restriction proj_;
		Proj2       proj2_;
		friend restriction;
		bind_transform_t(restriction proj, Proj2 proj2) : proj_{std::move(proj)}, proj2_{std::move(proj2)} {}

	 public:
		BOOST_MULTI_HD constexpr auto operator()(restriction::index idx) const noexcept { return proj2_(proj_[idx]); }
	};

 public:
	/// yields a restriction transformed in the leading dimension
	template<class Proj2, dimensionality_type One = 1 /*workaround for MSVC*/>
	BOOST_MULTI_HD auto transformed(Proj2 proj2) const -> restriction<1, bind_transform_t<Proj2>> {
		return bind_transform_t<Proj2>{*this, proj2} ^ multi::extents_t<One>({extent()});
	}

 private:
	template<class Fun, class... Args>
	static BOOST_MULTI_HD constexpr auto apply_(Fun&& fun, Args&&... args) {  // NOLINT(readability-identifier-naming) TODO(correaa) move this to ::detail
		using std::apply;
		return apply(std::forward<Fun>(fun), std::forward<Args>(args)...);
	}

	template<class Cursor, dimensionality_type DD = D>
	class cursor_t {
		Proj const* Pproj_;
		Cursor      cur_;
		friend class restriction;
		explicit constexpr cursor_t(Proj const* Pproj, Cursor cur) : Pproj_{Pproj}, cur_{cur} {}

	 public:
		using difference_type = restriction::difference_type;

		BOOST_MULTI_HD constexpr auto operator[](difference_type n) const -> decltype(auto) {
			if constexpr(DD != 1) {
				auto cur = cur_[n];
				return cursor_t<decltype(cur), DD - 1>{Pproj_, cur};
			} else {
				return apply_(*Pproj_, cur_[n]);
			}
		}
	};

 public:
	/// returns a cursor pointing to the top corner element of the array (a cursors is a lightweight representation of the array that drops the extents with pointer semantics)
	auto home() const {
		auto cur = extents().home();
		return cursor_t<decltype(cur), D>{&proj_, cur};
	}

	/// Random-access iterator in the leading dimension, in general they dereference to a restriction array of lower dimension or, for `D == 1`, to an element value (`T`) lazily generated.
	class iterator {
		using function_ptr = std::decay_t<decltype(&std::declval<Proj const&>())>;

		typename extents_t<D>::iterator it_;
		function_ptr                    Pproj_;
		// Proj const* Pproj_;

		iterator(typename extents_t<D>::iterator iter, function_ptr Pproj) : it_{iter}, Pproj_{Pproj} {}

		friend restriction;

	 public:
		constexpr iterator() = default;  // cppcheck-suppress uninitMemberVar ; partially formed

		iterator(iterator const& other) = default;
		iterator(iterator&&) noexcept   = default;

		auto operator=(iterator&&) noexcept -> iterator& = default;
		auto operator=(iterator const&) -> iterator&     = default;

		~iterator() = default;

	 private:
		using element_ref_  = typename detail::invoke_result_from_tuple<Proj, typename extents_t<D>::element>::type;
		using element_type_ = std::decay_t<element_ref_>;

	 public:
		/// Tag indicating the execution space (for compatibility with Thrust)
		using system = typename multi::detail::function_system<Proj>::type;

		/// A signed integer type for index arithmetic
		using difference_type = std::ptrdiff_t;
		/// An array of lower dimension or element type (e.g. `T`)
		using value_type      = std::conditional_t<
				 (D != 1),
				 restriction<D - 1, detail::bind_front_t<Proj>>,  // TODO(correaa) or multi::array
				 element_type_>;

		/// `void` (no pointer behind since it generates output)
		using pointer = void;

		/// `void` (the iterator generates output)
		using reference = std::conditional_t<(D != 1), restriction<D - 1, detail::bind_front_t<Proj>>, element_type_>;

		/// Iterator category (random access output)
		using iterator_category = std::random_access_iterator_tag;

		constexpr auto operator++() -> auto& {
			++it_;
			return *this;
		}
		constexpr auto operator--() -> auto& {
			--it_;
			return *this;
		}

		/// Addition assignment operator
		constexpr auto operator+=(difference_type diff) -> auto& {
			it_ += diff;
			return *this;
		}

		/// Substraction assignment operator
		constexpr auto operator-=(difference_type diff) -> auto& {
			it_ -= diff;
			return *this;
		}

		constexpr auto operator++(int) -> iterator {
			iterator ret{*this};
			++(*this);
			return ret;
		}
		constexpr auto operator--(int) -> iterator {
			iterator ret{*this};
			--(*this);
			return ret;
		}

		friend constexpr auto operator-(iterator const& self, iterator const& other) { return self.it_ - other.it_; }
		template<class = void>
		friend BOOST_MULTI_HD constexpr auto operator+(iterator const& self, difference_type n) {
			iterator ret{self};  // mull-ignore: cxx_init_const
			return ret += n;
		}
		friend constexpr auto operator-(iterator const& self, difference_type n) {
			iterator ret{self};
			return ret -= n;
		}

		friend constexpr auto operator+(difference_type n, iterator const& self) { return self + n; }

		friend constexpr auto operator==(iterator const& self, iterator const& other) noexcept -> bool { return self.it_ == other.it_; }
		friend constexpr auto operator!=(iterator const& self, iterator const& other) noexcept -> bool { return self.it_ != other.it_; }

		friend auto operator<=(iterator const& self, iterator const& other) noexcept -> bool { return self.it_ <= other.it_; }
		friend auto operator<(iterator const& self, iterator const& other) noexcept -> bool { return self.it_ < other.it_; }
		friend auto operator>(iterator const& self, iterator const& other) noexcept -> bool { return self.it_ > other.it_; }
		friend auto operator>=(iterator const& self, iterator const& other) noexcept -> bool { return self.it_ >= other.it_; }

		BOOST_MULTI_HD constexpr auto operator*() const -> reference {
			if constexpr(D != 1) {
				using std::get;
				// auto ll = [idx = get<0>(*it_), proj = proj_](auto... rest) { return proj(idx, rest...); };
				return restriction<D - 1, detail::bind_front_t<Proj>>(extents_t<D - 1>((*it_).tail()), detail::bind_front_t<Proj>{get<0>(*it_), *Pproj_});  // NOLINT(readability-redundant-parentheses) bug in clang-tidy trunk (May 2026)
			} else {
				using std::get;
#ifdef __CUDA_ARCH__
				return (const_cast<decltype(*Pproj_)&>(*Pproj_))(get<0>(*it_));  // workaround for __nv_dl_wrapper_t<__nv_dl_tag<int (*)(), main, 4>, int, int> >
#else
				return (*Pproj_)(get<0>(*it_));  // NOLINT(readability-redundant-parentheses) // bug in clang-tidy, can't be -> for CUDA
#endif
			}
		}

		BOOST_MULTI_HD constexpr auto operator[](difference_type diff) const -> reference { return *((*this) + diff); }  // TODO(correaa) use ra_iterator_facade
	};

	/// returns an iterator to the beginning in the leading dimension
	constexpr auto begin() const { return iterator{xs_.begin(), &proj_}; }
	/// returns an iterator to the end in the leading dimension
	constexpr auto end() const { return iterator{xs_.end(), &proj_}; }

	/// returns the size of the restriction array in the leading dimension.
	constexpr auto size() const noexcept { return xs_.size(); }

	/// returns a tuple decribing sizes across all the dimensions of the restriction array.
	constexpr auto sizes() const { return xs_.sizes(); }

	/// returns the extent of valid indices in the leading dimension.
	[[nodiscard]] constexpr auto extent() const { return xs_.extent(); }

	[[deprecated("use extent")]] constexpr auto extension() const { return xs_.extent(); }

	/// returns the index extents (structured cartesian product of index ranges) for all dimensions as an extents_type (`extents_t<D>`), a tuple of `D` index_extension values each encoding the extent in each dimension. The result can be passed directly to array constructors or compared for shape equality.
	[[nodiscard]] constexpr auto extents() const { return xs_; }

	/// Front subarray restriction in the leading dimensions or, in `D == 1`, the front element
	constexpr auto front() const { return *begin(); }
	/// Back subarray restriction in the leading dimensions or, in `D == 1`, the back element
	constexpr auto back() const { return *(begin() + (size() - 1)); }

 private:
	using elements_t = detail::restriction_elements_t<D, Proj>;

 public:
	/// yields a random‐access output range with all the elements of the array
	constexpr auto elements() const { return elements_t{xs_.elements(), proj_}; }

	/// returns the total number of elements in the array
	constexpr auto num_elements() const { return xs_.num_elements(); }
};

#ifdef __cpp_deduction_guides
/// deduces the dimensionality of the restriction from the extents and the restricted function
template<dimensionality_type D, typename Fun>
restriction(multi::extents_t<D>, Fun) -> restriction<D, Fun>;

template<typename Fun> restriction(extents_t<0>, Fun) -> restriction<0, Fun>;
template<typename Fun> restriction(extents_t<1>, Fun) -> restriction<1, Fun>;
/// Specialization for concrete extents arguments to allow a terse syntax, e.g. `auto r = restriction({3, 4}, fun);`.
template<typename Fun> restriction(extents_t<2>, Fun) -> restriction<2, Fun>;
template<typename Fun> restriction(extents_t<3>, Fun) -> restriction<3, Fun>;
template<typename Fun> restriction(extents_t<4>, Fun) -> restriction<4, Fun>;
template<typename Fun> restriction(extents_t<5>, Fun) -> restriction<5, Fun>;
template<typename Fun> restriction(extents_t<6>, Fun) -> restriction<6, Fun>;
#endif

template<dimensionality_type D, typename F>
auto restricted(F&& fun, extents_t<D> const& ext) {  // nvc++ has 'restrict' reserved
	return restriction<D, std::decay_t<F>>(ext, std::forward<F>(fun));
}

/// creates a restriction of `D` dimensions of a function that takes `D` arguments given an extents of Cartesian indices.
template<class F, dimensionality_type D>
BOOST_MULTI_HD constexpr auto operator^(F fun, extents_t<D> const& exts) {
	return restriction<D, F>(exts, std::move(fun));
}

/// Experimental namespace for the library
namespace experimental {
/// Binary operators that combines a D-dimensional extents object with a D-ary regularly invocable function object (generally a lambda)
template<class F, dimensionality_type D>
BOOST_MULTI_HD constexpr auto operator->*(extents_t<D> const& exts, F fun) {
	return restriction<D, F>(exts, std::move(fun));
}
}  // end namespace experimental

}  // end namespace boost::multi

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L) && !defined(_MSC_VER)
namespace std::ranges {  // NOLINT(cert-dcl58-cpp) to enable borrowed, nvcc needs namespace
template<class Fun, ::boost::multi::dimensionality_type D>
[[maybe_unused]] constexpr bool enable_borrowed_range<::boost::multi::restriction<D, Fun>> = true;  // NOLINT(misc-definitions-in-headers)
}  // end namespace std::ranges
#endif

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_RESTRICTION_HPP
