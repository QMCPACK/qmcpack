// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A reference implementation of a *hierarchical* `for_each` over a segmented range,
// in the spirit of Matt Austern's "Segmented Iterators and Hierarchical Algorithms" (1998).
//
// A segmented iterator denotes a position in a range that is internally a sequence of
// contiguous *segments* (e.g. a deque's blocks, a partitioned/distributed container, or
// `multi`'s row-major storage seen as rows of contiguous elements). Instead of walking the
// flattened range element-by-element, a hierarchical algorithm walks segment-by-segment
// (outer loop) and element-by-element within each segment (inner loop), which is both faster
// and exposes the contiguous inner runs to vectorization / specialized inner algorithms.
//
// To make a type `It` usable here, specialize `boost::multi::segmented_iterator_traits<It>`
// and provide the members documented on the primary template below.

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#include <algorithm>    // for std::for_each
#include <functional>   // IWYU pragma: keep  // for std::ref
#include <iterator>     // for std::next
#include <type_traits>  // for std::false_type, std::true_type
#include <utility>      // IWYU pragma: keep

namespace boost::multi {

// Customization point. The PRIMARY template marks a type as NOT segmented, so non-segmented
// (ordinary) iterators bottom out on a flat `std::for_each`. Specialize it for a segmented
// iterator type `It` and provide:
//
//   using is_segmented_iterator = std::true_type;
//   using iterator         = It;               // the segmented iterator itself
//   using segment_iterator = ...;              // iterates over the segments (the outer level)
//   using local_iterator   = ...;              // iterates the elements inside one segment (inner level)
//
//   static segment_iterator segment(iterator it);          // the segment that contains `it`
//   static local_iterator   local  (iterator it);          // position of `it` within its segment
//   static local_iterator   begin  (segment_iterator s);   // first element of segment `s`
//   static local_iterator   end    (segment_iterator s);   // one-past-the-last element of segment `s`
//   static iterator         compose(segment_iterator, local_iterator);  // (optional) inverse of segment()/local()
//
// `local_iterator` may itself be a segmented iterator: the algorithm recurses, so arbitrarily
// nested segmentation is handled. Contract: for a valid range [first, last), the segments
// `segment(first) .. segment(last)` and the local positions must be consistent so that the
// accesses `begin(segment(last))` and `local(last)` below are valid (the usual segmented-range
// invariant; e.g. an end position sits at `local == begin` of its segment, making the trailing
// partial segment empty).

namespace {
template<class Iterator>
struct segmented_iterator_traits {
	using is_segmented_iterator = std::false_type;
};
}  // namespace

namespace detail {
namespace {
// Applies `fun` (threaded BY REFERENCE) to every element of [first, last). If `It` is itself
// segmented, recurse with the hierarchical decomposition; otherwise bottom out on a flat
// `std::for_each`. Threading by reference (via std::ref at the leaf) lets a stateful, possibly
// move-only / non-assignable functor (e.g. a capturing lambda) accumulate across all segments.
template<class It, class UnaryFunction>
[[maybe_unused]] constexpr void for_each_segmented_apply(It first, It last, UnaryFunction& fun) {
	if constexpr(segmented_iterator_traits<It>::is_segmented_iterator::value) {
		using traits = segmented_iterator_traits<It>;

		auto const sfirst = traits::segment(first);
		auto const slast  = traits::segment(last);

		if(sfirst == slast) {  // the whole range lives in a single segment
			for_each_segmented_apply(traits::local(first), traits::local(last), fun);
			return;
		}

		// first segment (possibly partial): from the local position to the end of that segment
		for_each_segmented_apply(traits::local(first), traits::end(sfirst), fun);

		// middle segments (complete): begin() .. end() of each
		for(auto seg = std::next(sfirst); seg != slast; ++seg) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch)
			for_each_segmented_apply(traits::begin(seg), traits::end(seg), fun);
		}

		// last segment (possibly partial): from its begin to the local position
		for_each_segmented_apply(traits::begin(slast), traits::local(last), fun);
	} else {
		std::for_each(first, last, std::ref(fun));  // flat leaf; std::ref keeps mutating the same `fun`
	}
}
}  // namespace
}  // end namespace detail

namespace {
// Hierarchical for_each over a *segmented* range [first, last). Mirrors the contract of
// std::for_each: applies `fun` to every element exactly once (segment-by-segment), and returns
// the (moved-through) function object.
template<class SegmentedIterator, class UnaryFunction>
[[maybe_unused]] constexpr auto for_each_segmented(SegmentedIterator first, SegmentedIterator last, UnaryFunction fun) -> UnaryFunction {
	static_assert(
		segmented_iterator_traits<SegmentedIterator>::is_segmented_iterator::value,
		"for_each_segmented requires a segmented iterator (specialize boost::multi::segmented_iterator_traits<It>)"
	);
	detail::for_each_segmented_apply(first, last, fun);
	return fun;
}
}  // end namespace

}  // end namespace boost::multi

auto main() -> int {  // NOLINT(bugprone-exception-escape)
	// Intentionally empty: supply your own segmented iterator together with a
	// `boost::multi::segmented_iterator_traits<>` specialization, then test it against
	// `boost::multi::for_each_segmented` (e.g. compare its result to a flat std::for_each).
	return boost::report_errors();
}
