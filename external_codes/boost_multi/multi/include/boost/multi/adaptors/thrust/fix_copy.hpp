// Copyright 2021-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_THRUST_FIX_COPY_HPP
#define BOOST_MULTI_ADAPTORS_THRUST_FIX_COPY_HPP
#pragma once

namespace boost::multi {

#if 0
template<class Q1, class L1, class Size, class Q2, class R2, class L2>
auto copy_n(
	boost::multi::elements_iterator_t<                  Q1*                                                , L1>   first, Size count,
	boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::hip::tag, R2>, L2> d_first
)-> boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::hip::tag, R2>, L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		if constexpr(L1::dimensionality == 1 and L2::dimensionality == 1) {
			if(first.layout().stride() == 1 and d_first.layout().stride() == 1) {
				auto s = hipMemcpy  (raw_pointer_cast(d_first.current()),                                                                 first.current(),                                                               sizeof(Q2)* static_cast<std::size_t>(count), hipMemcpyHostToDevice); assert( s == hipSuccess );
			} else {
				auto s = hipMemcpy2D(raw_pointer_cast(d_first.current()), static_cast<std::size_t>(d_first.layout().stride())*sizeof(Q2), first.current(), static_cast<std::size_t>(first.layout().stride())*sizeof(Q2), sizeof(Q2), static_cast<std::size_t>(count), hipMemcpyHostToDevice); assert( s == hipSuccess );
			}
			return d_first + count;
		} else if constexpr(L1::dimensionality == 2 and L1::dimensionality == 2) {
			if(std::get<1>(first.layout().strides()) == 1 and std::get<1>(d_first.layout().strides()) == 1 and count%std::get<1>(first.layout().sizes()) == 0) {
				{auto s = hipMemcpy2D(raw_pointer_cast(d_first.current()), static_cast<std::size_t>(d_first.layout().stride())*sizeof(Q2), first.current(), static_cast<std::size_t>(first.layout().stride())*sizeof(Q2), static_cast<std::size_t>(std::get<1>(first.layout().sizes()))*sizeof(Q2), static_cast<std::size_t>(count/std::get<1>(first.layout().sizes())), hipMemcpyHostToDevice); assert( s == cudaSuccess );}
				return d_first + count;
			}  // else fallthrough
		}
		{auto r = hipHostRegister(
			const_cast<void*>(static_cast<void const*>(first.base())),
			static_cast<std::size_t>                  (first.layout().hull_size()*sizeof(Q1)),
			hipHostRegisterPortable
		); assert( r == hipSuccess );}
		auto ret = ::thrust::copy_n(
			::thrust::hip::par,
			first, count, d_first
		);
		{auto r = hipHostUnregister(
			const_cast<void*>(static_cast<void const*>(first.base()))
		); assert( r == hipSuccess );}
		return ret;
	} else {
		return ::thrust::copy_n(first, count, d_first);
	}
	return d_first + count;
}

template<class Q1, class R1, class L1, class Size, class Q2, class L2>
auto copy_n(
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::hip::tag, R1, ::thrust::use_default>, L1>   first, Size count,
	boost::multi::elements_iterator_t<                  Q2*                                                    , L2> d_first
)-> boost::multi::elements_iterator_t<                  Q2*                                                    , L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		if constexpr(L1::dimensionality == 1 and L2::dimensionality == 1) {
			if(first.layout().stride() == 1 and d_first.layout().stride() == 1) {
				auto s = hipMemcpy  (                 d_first.current() ,                                                                 raw_pointer_cast(first.current()),                                                               sizeof(Q2)* static_cast<std::size_t>(count), hipMemcpyDeviceToHost); assert( s == cudaSuccess );
			} else {
				auto s = hipMemcpy2D(                 d_first.current() , static_cast<std::size_t>(d_first.layout().stride())*sizeof(Q2), raw_pointer_cast(first.current()), static_cast<std::size_t>(first.layout().stride())*sizeof(Q2), sizeof(Q2), static_cast<std::size_t>(count), hipMemcpyDeviceToHost); assert( s == cudaSuccess );
			}
			return d_first + count;
		} else if constexpr(L1::dimensionality == 2 and L1::dimensionality == 2) {
			if(std::get<1>(first.layout().strides()) == 1 and std::get<1>(d_first.layout().strides()) == 1 and count%std::get<1>(first.layout().sizes()) == 0) {
				auto s = hipMemcpy2D(                 d_first.current() , static_cast<std::size_t>(d_first.layout().stride())*sizeof(Q2), raw_pointer_cast(first.current()), static_cast<std::size_t>(first.layout().stride())*sizeof(Q2), static_cast<std::size_t>(std::get<1>(first.layout().sizes()))*sizeof(Q2), static_cast<std::size_t>(count/std::get<1>(first.layout().sizes())), hipMemcpyDeviceToHost); assert( s == cudaSuccess );
				return d_first + count;
			}
		}
		{auto r = hipHostRegister(
			const_cast<void*>(static_cast<void const*>(d_first.base())),
			static_cast<std::size_t>                  (d_first.layout().hull_size()*sizeof(Q1)),
			hipHostRegisterPortable
		); assert( r == hipSuccess );}
		auto ret = ::thrust::copy_n(
			::thrust::hip::par,
			first, count, d_first
		);
		{auto r = hipHostUnregister(
			const_cast<void*>(static_cast<void const*>(d_first.base()))
		); assert( r == hipSuccess );}
		return ret;
	} else {
		return ::thrust::copy_n(first, count, d_first);
	}
	return d_first + count;
}

template<class Q1, class L1, class Size, class Q2, class R2, class L2>
auto uninitialized_copy_n(
	boost::multi::elements_iterator_t<                  Q1*                                                    , L1>   first, Size count,
	boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::hip::tag, R2, ::thrust::use_default>, L2> d_first
)-> boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::hip::tag, R2, ::thrust::use_default>, L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		return boost::multi::copy_n(first, count, d_first);
	} else {
		return ::thrust::uninitialized_copy_n(first, count, d_first);
	}
}

template<class Q1, class R1, class L1, class Size, class Q2, class L2>
auto uninitialized_copy_n(
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::hip::tag, R1, ::thrust::use_default>, L1>   first, Size count,
	boost::multi::elements_iterator_t<                  Q2*                                                    , L2> d_first
)-> boost::multi::elements_iterator_t<                  Q2*                                                    , L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		return boost::multi::copy_n(first, count, d_first);
	} else {
		return ::thrust::uninitialized_copy_n(first, count, d_first);
	}
}

template<class Q1, class L1, class Q2, class R2, class L2>
auto copy(
	boost::multi::elements_iterator_t<                  Q1*                             , L1>   first,
	boost::multi::elements_iterator_t<                  Q1*                             , L1>   last ,
	boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::hip::tag, R2>, L2> d_first
)-> boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::hip::tag, R2>, L2> {
	return boost::multi::copy_n(first, last - first, d_first);
}

template<class Q1, class R1, class L1, class Q2, class L2>
auto copy(
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::hip::tag, R1>, L1>   first,
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::hip::tag, R1>, L1>   last ,
	boost::multi::elements_iterator_t<                  Q2*                             , L2> d_first
)-> boost::multi::elements_iterator_t<                  Q2*                             , L2> {
	return boost::multi::copy_n(first, last - first, d_first);
}

template<class Q1, class L1, class Q2, class R2, class L2>
auto uninitialized_copy(
	boost::multi::elements_iterator_t<                  Q1*                             , L1>   first,
	boost::multi::elements_iterator_t<                  Q1*                             , L1>   last ,
	boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::hip::tag, R2>, L2> d_first
)-> boost::multi::elements_iterator_t<::thrust::pointer<Q2, ::thrust::hip::tag, R2>, L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		return boost::multi::copy(first, last, d_first);
	} else {
		return ::thrust::uninitialized_copy(first, last, d_first);
	}
}

template<class Q1, class R1, class L1, class Q2, class L2>
auto uninitialized_copy(
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::hip::tag, R1>, L1>   first,
	boost::multi::elements_iterator_t<::thrust::pointer<Q1, ::thrust::hip::tag, R1>, L1>   last ,
	boost::multi::elements_iterator_t<                  Q2*                             , L2> d_first
)-> boost::multi::elements_iterator_t<                  Q2*                             , L2> {
	if constexpr(std::is_trivially_assignable<Q2&, Q1&>{}) {
		return boost::multi::copy(first, last, d_first);
	} else {
		return ::thrust::uninitialized_copy(first, last, d_first);
	}
}

#endif

}  // end namespace boost::multi

#endif  // BOOST_MULTI_ADAPTORS_THRUST_FIX_COPY_HPP
