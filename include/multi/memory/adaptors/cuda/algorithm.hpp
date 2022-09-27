#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXXX $CXXFLAGS $0 -o $0x -lcudart -lboost_unit_test_framework -lboost_timer&&$0x&&rm $0x;exit
#endif
#ifndef BOOST_MULTI_MEMORY_ADAPTORS_CUDA_ALGORITHM_HPP
#define BOOST_MULTI_MEMORY_ADAPTORS_CUDA_ALGORITHM_HPP

#include          "../cuda/cstring.hpp"
#include "../../../array_ref.hpp"
#include "../../../config/MARK.hpp"
#include "../../../detail/adl.hpp"

#include <thrust/device_ptr.h>
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include "../cuda/error.hpp"

#include<complex> //TODO remove, handle complex_fix instead?

namespace boost {namespace multi {
namespace memory {namespace cuda {

template<class U, class T, typename Size, typename = std::enable_if_t<std::is_trivially_assignable<T&, U>{}>>
ptr<T> copy_n(ptr<U> first, Size count, ptr<T> result) {
	return memcpy(result, first, count*sizeof(T)), result + count;
}

template<class U, class T, typename Size, typename = std::enable_if_t<std::is_trivially_assignable<T&, U>{}>>
ptr<T> copy_n(U* first, Size count, ptr<T> result) {return memcpy(result, first, count*sizeof(T)), result + count;}

template<class U, class T, typename Size, typename = std::enable_if_t<std::is_trivially_assignable<T&, U>{}>>
T*     copy_n(ptr<U> first, Size count, T* result) {return memcpy(result, first, count*sizeof(T)), result + count;}

template<class U, class T, typename Size, typename = std::enable_if_t<std::is_trivially_assignable<T&, U>{}>>
ptr<T> uninitialized_copy_n(ptr<U> first, Size count, ptr<T> result) {return memcpy(result, first, count*sizeof(T)), result + count;}

template<class U, class T, typename Size, typename = std::enable_if_t<std::is_trivially_assignable<T&, U>{}>>
ptr<T> uninitialized_copy_n(U* first, Size count, ptr<T> result) {return memcpy(result, first, count*sizeof(T)), result + count;}

template<class U, class T, typename Size, typename = std::enable_if_t<std::is_trivially_assignable<T&, U>{}>>
T* uninitialized_copy_n(ptr<U> first, Size count, T* result) {return memcpy(result, first, count*sizeof(T)), result + count;}

template<class PtrU, class T>
auto copy(PtrU first, PtrU last, ptr<T> result){
	return copy_n(first, std::distance(first, last), result);
}

template<class U, class T>
auto copy(ptr<U> first, ptr<U> last, ptr<T> result){
	return copy_n(first, std::distance(first, last), result);
}

//->decltype(copy_n(first, std::distance(first, last), result)){
//	return copy_n(first, std::distance(first, last), result);}


template<class T>
auto fill(memory::cuda::ptr<T> first, memory::cuda::ptr<T> last, T const& value)
->decltype(fill_n(first, std::distance(first, last), value)){
	return fill_n(first, std::distance(first, last), value);}

template<class U, class T, typename Size, typename = std::enable_if_t<std::is_trivially_assignable<T&, U>{}>>
memory::cuda::ptr<T> fill_n(ptr<T> const first, Size count, U const& value){
	if(std::find_if((char const*)(&value), (char const*)(&value) + sizeof(value), [](char c){return c!=0;}) == (char const*)(&value) + sizeof(value)){
//	if(std::find(reinterpret_cast<char const*>(&value), reinterpret_cast<char const*>(&value) + sizeof(value), true) == reinterpret_cast<char const*>(&value) + sizeof(value)){
//	if(value == 0.){
		cuda::memset(first, 0, count*sizeof(T));
	}
	else if(count--) for(ptr<T> new_first = adl_copy_n(&value, 1, first); count;){
		auto n = std::min(Size(std::distance(first, new_first)), count);
		new_first = copy_n(first, n, new_first);
		count -= n;
	}
	return first + count;
}

template<class U, class T, class TT, typename Size>//, typename = std::enable_if_t<std::is_trivially_assignable<T&, U>{}>>
memory::cuda::ptr<T> fill_n(array_iterator<T, 1, ptr<TT>> const first, Size count, U const& value){
	if(count--) 
		for(ptr<T> new_first = adl_copy_n(&value, 1, first); count;){
			auto n = std::min(Size(std::distance(first, new_first)), count);
			new_first = copy_n(first, n, new_first);
			count -= n;
		}
	return first + count;
}

// TODO fix traits std::complex instead
template<class T1, typename Size, class T2, std::enable_if_t<std::is_same<std::decay_t<T1>, T2>{},int> =0>
ptr<std::complex<T2>> copy_n(ptr<T1> first, Size count, ptr<std::complex<T2>> result){
	fill_n(result, count, std::complex<T2>{0});
	copy_n(
		multi::array_iterator<std::decay_t<T1>, 1, ptr<T1>>{first, 1}, count,
		multi::array_iterator<T2, 1, ptr<T2>>{reinterpret_pointer_cast<T2>(result), 2}
	);
	return result + count;
}

template<class T1, class TT1, typename Size, class T2, std::enable_if_t<std::is_same<std::decay_t<T1>, T2>{},int> =0>
auto copy_n(iterator<T1, 1, ptr<TT1>> first, Size count, iterator<std::complex<T2>, 1, ptr<std::complex<T2>>> result){
	MULTI_MARK_SCOPE("cuda copy_n 1D complex destination");
	if(stride(first) == 1 and stride(result)==1) copy_n(base(first), count, base(result));
	else assert(0);
	return result + count;
}

//template<class T, class Size, class V, std::enable_if_t<std::is_trivially_default_constructible<T>{}, int> =0>
//auto uninitialized_fill_n(cuda::ptr<T> first, Size n, V const& v){return fill_n(first, n, v);}

template<class TT, class T, class Size, std::enable_if_t<std::is_trivially_default_constructible<typename std::iterator_traits<ptr<TT>>::value_type>{}, int> =0>
auto uninitialized_fill_n(ptr<TT> first, Size n, T const& t)
{
	return fill_n(first, n, t);
}


template<class It, class Size, class T = typename std::iterator_traits<It>::value_type, std::enable_if_t<std::is_trivially_constructible<T>{}, int> = 0>
auto uninitialized_value_construct_n(It first, Size n){
	return uninitialized_fill_n(first, n, T());
}

template<class It, class Size, class T, std::enable_if_t<std::is_trivially_default_constructible<typename std::iterator_traits<ptr<T>>::value_type>{}, int> = 0>
auto uninitialized_copy_n(It first, Size n, ptr<T> d_first)
->decltype(copy_n(first, n, d_first)){
	return copy_n(first, n, d_first);}

//template<class Alloc, class It, class Size>
//auto alloc_uninitialized_value_construct_n(Alloc&, It first, Size n){
//	return uninitialized_value_construct_n(first, n);
//}

#if 1

#if 1

template<int A, int B, int C, int D> void what() = delete;

template<class T1, class Q1, class Size, class T2, class Q2>
array_iterator<T2, 1, ptr<Q2>> copy_n(
	array_iterator<T1, 1, ptr<Q1>> first_ , Size count, 
	array_iterator<T2, 1, ptr<Q2>> result_
){
	MULTI_MARK_SCOPE("cuda copy_n 1D");
	array_iterator<T1, 1, ::thrust::device_ptr<Q1>> first ; std::memcpy((void*)&first , (void const*)&first_ , sizeof(first_));
	array_iterator<T2, 1, ::thrust::device_ptr<Q2>> result; std::memcpy((void*)&result, (void const*)&result_, sizeof(first_));
	static_assert( sizeof(first ) == sizeof(first_ ) );
	static_assert( sizeof(result) == sizeof(result_) );
	::thrust::for_each(
		::thrust::make_counting_iterator(0L),
		::thrust::make_counting_iterator(count),
		[first, result, x = multi::extensions_t<1>(count)] __device__ (auto n){ // requires --extended-lambda nvcc flag
		//  std::tuple<index>
			auto const i = x.from_linear(n);
			using boost::multi::detail::get;
			result[get<0>(i)] = T2(first[get<0>(i)]);
		}
	);
	return result_ + count;
}

template<class T1, class Q1, class Size, class T2, class Q2>
array_iterator<T2, 2, ptr<Q2>>
copy_n(
	array_iterator<T1, 2, ptr<Q1>> first_ , Size count,
	array_iterator<T2, 2, ptr<Q2>> result_
) {
	MULTI_MARK_SCOPE("cuda copy_n 2D");

	array_iterator<T1, 2, ::thrust::device_ptr<Q1>> first ; std::memcpy((void*)&first , (void const*)&first_ , sizeof(first_));
	array_iterator<T2, 2, ::thrust::device_ptr<Q2>> result; std::memcpy((void*)&result, (void const*)&result_, sizeof(first_));

	static_assert( sizeof(first ) == sizeof(first_ ) );
	static_assert( sizeof(result) == sizeof(result_) );

	assert(first->extensions() == result->extensions());
	::thrust::for_each(
		::thrust::make_counting_iterator(0L),
		::thrust::make_counting_iterator(count*first->num_elements()),
		[first, count, result, x = first->extensions()] __device__ (auto n){
			auto const ij = (count*x).from_linear(n);
			result[std::get<0>(ij)][std::get<1>(ij)] = T2(first[std::get<0>(ij)][std::get<1>(ij)]);
		}
	);

	return result_ + count;
}

template<class T1, class Q1, class Size, class T2, class Q2>
array_iterator<T2, 3, ptr<Q2>>
copy_n(
	array_iterator<T1, 3, ptr<Q1>> first_ , Size count, 
	array_iterator<T2, 3, ptr<Q2>> result_
) {
	MULTI_MARK_SCOPE("cuda copy_n 3D");
	array_iterator<T1, 3, ::thrust::device_ptr<Q1>> first ; std::memcpy((void*)&first , (void const*)&first_ , sizeof(first_));
	array_iterator<T2, 3, ::thrust::device_ptr<Q2>> result; std::memcpy((void*)&result, (void const*)&result_, sizeof(first_));
	static_assert( sizeof(first ) == sizeof(first_ ) );
	static_assert( sizeof(result) == sizeof(result_) );
	assert(first->extensions() == result->extensions());
	::thrust::for_each(
		::thrust::make_counting_iterator(0L),
		::thrust::make_counting_iterator(count*first->num_elements()),
		[first, count, result, x = first->extensions()] __device__ (auto n){
			auto const ijk = (count*x).from_linear(n);
			result.apply(ijk) = T2(first.apply(ijk));
		}
	);
	return result_ + count;
}

template<class T1, class Q1, class Size, class T2, class Q2>
array_iterator<T2, 4, ptr<Q2>>
copy_n(
	array_iterator<T1, 4, ptr<Q1>> first_ , Size count, 
	array_iterator<T2, 4, ptr<Q2>> result_
) {
	MULTI_MARK_SCOPE("cuda copy_n 4D");
	array_iterator<T1, 4, ::thrust::device_ptr<Q1>> first ; std::memcpy((void*)&first , (void const*)&first_ , sizeof(first_));
	array_iterator<T2, 4, ::thrust::device_ptr<Q2>> result; std::memcpy((void*)&result, (void const*)&result_, sizeof(first_));
	static_assert( sizeof(first ) == sizeof(first_ ) );
	static_assert( sizeof(result) == sizeof(result_) );
	assert(first->extensions() == result->extensions());
	::thrust::for_each(
		::thrust::make_counting_iterator(0L),
		::thrust::make_counting_iterator(count*first->num_elements()),
		[first, count, result, x = first->extensions()] __device__ (auto n){
			auto const ijk = (count*x).from_linear(n);
			result.apply(ijk) = T2(first.apply(ijk));
		}
	);
	return result_ + count;
}

template<class T1, class Q1, class Size, class T2, class Q2>
array_iterator<T2, 5, ptr<Q2>> 
copy_n(
	array_iterator<T1, 5, ptr<Q1>> first_ , Size count, 
	array_iterator<T2, 5, ptr<Q2>> result_
) {
	MULTI_MARK_SCOPE("cuda copy_n 5D");
	array_iterator<T1, 5, ::thrust::device_ptr<Q1>> first ; std::memcpy((void*)&first , (void const*)&first_ , sizeof(first_));
	array_iterator<T2, 5, ::thrust::device_ptr<Q2>> result; std::memcpy((void*)&result, (void const*)&result_, sizeof(first_));
	static_assert( sizeof(first ) == sizeof(first_ ) );
	static_assert( sizeof(result) == sizeof(result_) );
	assert(first->extensions() == result->extensions());
	::thrust::for_each(
		::thrust::make_counting_iterator(0L),
		::thrust::make_counting_iterator(count*first->num_elements()),
		[first, count, result, x = first->extensions()] __device__ (auto n){
			auto const ijk = (count*x).from_linear(n);
			result.apply(ijk) = T2(first.apply(ijk));
		}
	);
	return result_ + count;
}


template<class T1, class Q1, class Size, class T2, class Q2, multi::dimensionality_type D, std::enable_if_t<(D<=5), int>* =0>
array_iterator<T2, D, ptr<Q2>> 
copy(
	array_iterator<T1, D, ptr<Q1>> first_ , array_iterator<T1, D, ptr<Q1>> last_, 
	array_iterator<T2, D, ptr<Q2>> result_
) {return copy_n(first_, last_ - first_, result_);}

#endif

template<class T1, class Q1, typename Size, class T2, class Q2, typename = std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}>>
auto copy_n(iterator<T1, 1, Q1*> first, Size count, iterator<T2, 1, ptr<Q2>> result)
->decltype(memcpy2D(result.base(), sizeof(T2)*stride(result), first.base(), sizeof(T1)*stride(first ), sizeof(T1), count), result + count){
	MULTI_MARK_SCOPE("cuda copy_n 1D cpu source gpu destination");
	return memcpy2D(base(result), sizeof(T2)*stride(result), base(first), sizeof(T1)*stride(first ), sizeof(T1), count), result + count;}

template<class T1, class Q1, typename Size, class T2, class Q2, typename = std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}>>
auto copy_n(iterator<T1, 1, ptr<Q1>> first, Size count, iterator<T2, 1, Q2*> result)
->decltype(memcpy2D(result.base(), sizeof(T2)*stride(result), first.base(), sizeof(T1)*stride(first), sizeof(T1), count), result+count){
	MULTI_MARK_SCOPE("cuda copy_n 1D gpu source cpu destination");
	return memcpy2D(base(result), sizeof(T2)*stride(result), base(first), sizeof(T1)*stride(first), sizeof(T1), count), result+count;}

//template<class T1, class Q1, typename Size, class T2, class Q2, typename = std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}>>
//auto copy_n(iterator<T1, 1, ptr<Q1>> first, Size count, iterator<T2, 1, ptr<Q2>> result)
//->decltype(memcpy2D(result.base(), sizeof(T2)*stride(result), first.base(), sizeof(T1)*stride(first), sizeof(T1), count), result + count){assert(0);
//	return memcpy2D(base(result), sizeof(T2)*stride(result), base(first), sizeof(T1)*stride(first), sizeof(T1), count), result + count;}

template<class T1, class... Q1, typename Size, class T2, class... Q2, typename = std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}>>
auto copy_n(iterator<T1, 1, managed::ptr<Q1...>> first, Size count, iterator<T2, 1, managed::ptr<Q2...>> result)
->decltype(memcpy2D(result.base(), sizeof(T2)*stride(result), first.base(), sizeof(T1)*stride(first), sizeof(T1), count), result + count){
	MULTI_MARK_SCOPE("cuda copy_n 1D gpu managed source gpu managed destination");
	return memcpy2D(base(result), sizeof(T2)*stride(result), base(first), sizeof(T1)*stride(first), sizeof(T1), count), result + count;}

template<class T1, class... Q1, class T2, class... Q2, typename = std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}>>
auto copy(iterator<T1, 1, managed::ptr<Q1...>> first, iterator<T1, 1, managed::ptr<Q1...>> last, iterator<T2, 1, managed::ptr<Q2...>> result)
->decltype(copy_n(first, last - first, result)){
	return copy_n(first, last - first, result);}

#endif

#define ENABLE_IF class=std::enable_if_t

template<class T1, class T1P, class Size, class T2, class T2P, 
//	std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}, int> =0
	ENABLE_IF<std::is_trivially_assignable<T2&, T1>{}>
//	typename = std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}>
>
auto copy_n(array_iterator<T1, 1, ptr<T1P>> first, Size count, array_iterator<T1, 1, ptr<T1P>> result)
//->decltype(memcpy2D(base(result), sizeof(T2)*stride(result), base(first), sizeof(T1)*stride(first), sizeof(T1), count), result + count)
{	return memcpy2D(base(result), sizeof(T2)*stride(result), base(first), sizeof(T1)*stride(first), sizeof(T1), count), result + count;}

template<class T1, class T1P, class T2, class T2P>
auto copy(array_iterator<T1, 1, cuda::ptr<T1P>> first, array_iterator<T1, 1, cuda::ptr<T1P>> last, array_iterator<T2, 1, cuda::ptr<T2P>> result)
//->decltype(cuda::copy_n(first, last - first, result))
{	return cuda::copy_n(first, last - first, result);}

//template<class T1, class P1, class T2, class P2>
//auto copy(array_iterator<T1, 1, cuda::ptr<P1 >> first, array_iterator<T1, 1, cuda::ptr<P1 >> last, array_iterator<T2, 1, cuda::ptr<P2>> result){
//	std::cout << boost::stacktrace::stacktrace() << std::endl;
//	assert(0);
//}


template<class T, class Size, std::enable_if_t<std::is_trivially_default_constructible<T>{}, int> =0>
auto uninitialized_value_construct_n(ptr<T> first, Size n){return uninitialized_fill_n(first, n, T{});}

namespace managed {

//template<class T1, class T1const, class... A1, class Size, class T2, class T2const, class... A2, typename = std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}>>
//auto copy_n(
//	array_iterator<T1, 1, managed::ptr<T1const, A1...>> first, Size n, 
//	array_iterator<T2, 1, managed::ptr<T2const, A2...>> result
//)
//->decltype(memcpy2D(base(result), sizeof(T2)*stride(result), base(first), sizeof(T1)*stride(first), sizeof(T1), n), result + n){
//	return memcpy2D(base(result), sizeof(T2)*stride(result), base(first), sizeof(T1)*stride(first), sizeof(T1), n), result + n;}

//template<class T1, class... A1, class Size, class T2, class... A2>//, std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}, int> =0>
//auto copy_n(
//	managed::ptr<T1, A1...> first, Size count, 
//	managed::ptr<T2, A2...> result
//)
//->decltype(cuda::copy_n(cuda::ptr<T1>(first), count, cuda::ptr<T2>(result))){
//{	return cuda::copy_n(cuda::ptr<T1>(first), count, cuda::ptr<T2>(result));}


template<class T1, class... A1, class Size, class T2, class... A2>//, std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}, int> =0>
auto copy_n(
	managed::ptr<T1, A1...> first, Size count, 
	managed::ptr<T2, A2...> result
)
//->decltype(cuda::copy_n(cuda::ptr<T1>(first), count, cuda::ptr<T2>(result))){
{	return cuda::copy_n(cuda::ptr<T1>(first), count, cuda::ptr<T2>(result));}

template<class T1, class T2>//, std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}, int> =0>
auto copy(
	managed::ptr<T1> first, managed::ptr<T1> last, 
	managed::ptr<T2> result
)
//->decltype(cuda::copy(first, last, result)){assert(0);
{	return cuda::copy(cuda::ptr<T1>(first), cuda::ptr<T1>(last), cuda::ptr<T2>(result)), result + (last - first);}


//template<class T1, typename Size, class T2, std::enable_if_t<std::is_same<std::decay_t<T1>, T2>{},int> =0>
inline
managed::ptr<std::complex<double>> copy_n(managed::ptr<double> /*first*/, std::size_t /*count*/, managed::ptr<std::complex<double>> result) {
	assert(0);
	return result;
}


template<class T1, class T1const, class... A1, class Size, class T2, class T2const, class... A2>
auto copy_n(
	array_iterator<T1, 1, managed::ptr<T1const, A1...>> first, Size count,
	array_iterator<T2, 1, managed::ptr<T2const, A2...>> d_first
){
	return cuda::copy_n(array_iterator<T1, 1, cuda::ptr<T1const, A1...>>(first), count, array_iterator<T2, 1, cuda::ptr<T2const, A2...>>(d_first)), d_first + count;
}

template<class T1, class T1const, class... A1, class Size, class T2, class T2const, class... A2>
auto copy_n(
	array_iterator<T1, 2, managed::ptr<T1const, A1...>> first, Size count,
	array_iterator<T2, 2, managed::ptr<T2const, A2...>> d_first
){	
	return cuda::copy_n(array_iterator<T1, 2, cuda::ptr<T1const, A1...>>(first), count, array_iterator<T2, 2, cuda::ptr<T2const, A2...>>(d_first)), d_first + count;
}

template<class T1, class T1const, class... A1, class Size, class T2, class T2const, class... A2>
auto copy_n(
	array_iterator<T1, 3, managed::ptr<T1const, A1...>> first, Size count,
	array_iterator<T2, 3, managed::ptr<T2const, A2...>> d_first
){	
	return cuda::copy_n(array_iterator<T1, 3, cuda::ptr<T1const, A1...>>(first), count, array_iterator<T2, 3, cuda::ptr<T2const, A2...>>(d_first)), d_first + count;
}

template<class T1, class T1const, class... A1, class Size, class T2, class T2const, class... A2>
auto copy_n(
	array_iterator<T1, 4, managed::ptr<T1const, A1...>> first, Size count,
	array_iterator<T2, 4, managed::ptr<T2const, A2...>> d_first
){	
	return cuda::copy_n(array_iterator<T1, 4, cuda::ptr<T1const, A1...>>(first), count, array_iterator<T2, 4, cuda::ptr<T2const, A2...>>(d_first)), d_first + count;
}

template<class T1, class T1const, class... A1, class Size, class T2, class T2const, class... A2>
auto copy_n(
	array_iterator<T1, 5, managed::ptr<T1const, A1...>> first, Size count,
	array_iterator<T2, 5, managed::ptr<T2const, A2...>> d_first
){	
	return cuda::copy_n(array_iterator<T1, 5, cuda::ptr<T1const, A1...>>(first), count, array_iterator<T2, 5, cuda::ptr<T2const, A2...>>(d_first)), d_first + count;
}


template<class T1, class T1const, class... A1, class T2, class T2const, class... A2>
auto copy(
	array_iterator<T1, 1, managed::ptr<T1const, A1...>> first,
	array_iterator<T1, 1, managed::ptr<T1const, A1...>> last,
	array_iterator<T2, 1, managed::ptr<T2const, A2...>> d_first
){
	return managed::copy_n(first, last - first, d_first);
}

template<class T1, class T1const, class... A1, class T2, class T2const, class... A2>
auto copy(
	array_iterator<T1, 2, managed::ptr<T1const, A1...>> first,
	array_iterator<T1, 2, managed::ptr<T1const, A1...>> last,
	array_iterator<T2, 2, managed::ptr<T2const, A2...>> d_first
){
	return managed::copy_n(first, last - first, d_first);
}

template<class T1, class T1const, class... A1, class T2, class T2const, class... A2>
auto copy(
	array_iterator<T1, 3, managed::ptr<T1const, A1...>> first,
	array_iterator<T1, 3, managed::ptr<T1const, A1...>> last,
	array_iterator<T2, 3, managed::ptr<T2const, A2...>> d_first
){
	return managed::copy_n(first, last - first, d_first);
}

template<class T1, class T1const, class... A1, class T2, class T2const, class... A2>
auto copy(
	array_iterator<T1, 4, managed::ptr<T1const, A1...>> first,
	array_iterator<T1, 4, managed::ptr<T1const, A1...>> last,
	array_iterator<T2, 4, managed::ptr<T2const, A2...>> d_first
){
	return managed::copy_n(first, last - first, d_first);
}

template<class T1, class T1const, class... A1, class T2, class T2const, class... A2>
auto copy(
	array_iterator<T1, 5, managed::ptr<T1const, A1...>> first,
	array_iterator<T1, 5, managed::ptr<T1const, A1...>> last,
	array_iterator<T2, 5, managed::ptr<T2const, A2...>> d_first
){
	return managed::copy_n(first, last - first, d_first);
}

template<class U, class T, class... As, typename Size>//, typename = std::enable_if_t<std::is_trivially_assignable<T&, U>{}>>
auto fill_n(managed::ptr<T, As...> first, Size count, U const& value){
	return fill_n(cuda::ptr<T>(first), count, value), first + count;
}

template<class T, class...As, class Size, class V, std::enable_if_t<std::is_trivially_constructible<T, V const&>{}, int> =0>
auto uninitialized_fill_n(managed::ptr<T, As...> first, Size n, V const& v){return fill_n(first, n, v);}

template<class Alloc, class Ptr, class ForwardIt, std::enable_if_t<std::is_trivially_copy_constructible<typename std::iterator_traits<ForwardIt>::value_type>{}, int> = 0>
auto alloc_uninitialized_copy(Alloc&, Ptr first, Ptr last, ForwardIt dest)
->decltype(cuda::copy(first, last, dest)){
	return cuda::copy(first, last, dest);}

template<class Alloc, class Ptr, class Size, class ForwardIt, std::enable_if_t<std::is_trivially_copyable<typename std::iterator_traits<ForwardIt>::value_type>{}, int> = 0>
auto alloc_uninitialized_copy(Alloc&, Ptr first, Size n, ForwardIt dest)
->decltype(cuda::copy_n(first, n, dest)){
	return cuda::copy_n(first, n, dest);}

template<class T, class TP, class Size, std::enable_if_t<std::is_trivially_default_constructible<T>{}, int> =0>
auto uninitialized_default_construct_n(cuda::managed::ptr<T, TP> first, Size n){return first + n;}

template<class T, class TP, class Size, std::enable_if_t<std::is_trivially_default_constructible<T>{}, int> =0>
auto uninitialized_default_construct_n(cuda::managed::ptr<std::complex<T>, TP> first, Size n){
	return first + n; // TODO remove?
}

template<class T, class TP, class Size, std::enable_if_t<std::is_trivially_destructible<T>{}, int> =0>
auto destroy_n(cuda::managed::ptr<T, TP> first, Size n){return first + n;}

}

}}

/*
template<class T1, typename Size, class T2, class Q2, typename = std::enable_if_t<std::is_trivially_assignable<T2&, T1>{}>>
array_iterator<T2, 1, memory::cuda::ptr<Q2>> copy_n(
	T1* first, Size count, 
	array_iterator<T2, 1, memory::cuda::ptr<Q2>> result
){
	std::cout << "count " << std::endl;
	return copy_n<array_iterator<std::decay_t<T1>, 1, T1*>>(first, count, result);
}*/

#if 0
template<class T1, class Q1, class T2, class Q2>
auto copy(
	array_iterator<T1, 1, memory::cuda::ptr<Q1>> f, array_iterator<T1, 1, memory::cuda::ptr<Q1>> l, 
	array_iterator<T2, 1, memory::cuda::ptr<Q2>> d
)
->decltype(copy_n(f, std::distance(f, l), d)){assert(stride(f)==stride(l));
	return copy_n(f, std::distance(f, l), d);}

template<class T1, class Q1, class T2, class Q2>
auto copy(
	array_iterator<T1, 1, Q1*> f, array_iterator<T1, 1, memory::cuda::ptr<Q1>> l, 
	array_iterator<T2, 1, memory::cuda::ptr<Q2>> d
)
->decltype(copy_n(f, std::distance(f, l), d)){assert(stride(f)==stride(l));
	return copy_n(f, std::distance(f, l), d);}
#endif

}}

#define FWD(x) std::forward<decltype(x)>(x)

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//#if not __INCLUDE_LEVEL__ // _TEST_MULTI_MEMORY_ADAPTORS_CUDA_ALGORITHM

//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi initializer_list"
//#define BOOST_TEST_DYN_LINK
//#include<boost/test/unit_test.hpp>

//#include "../../../array.hpp"
//#include "../../../adaptors/cuda.hpp"

//#include "../cuda/allocator.hpp"
//#include<boost/timer/timer.hpp>
//#include<numeric>

//namespace multi = boost::multi;
//namespace cuda = multi::memory::cuda;

//BOOST_AUTO_TEST_CASE(copy_1d){
//	auto const A_cpu = []{
//		multi::array<double, 1> ret(10);
//		std::generate(
//			ret.data_elements(), ret.data_elements() + ret.num_elements(), std::rand
//		);
//		return ret;
//	}();
//	std::cout<<"memory size "<< A_cpu.num_elements()*sizeof(decltype(A_cpu)::element)/1e6 <<" MB\n";
//	{
//		multi::array<double, 1> B_cpu(size(A_cpu));
//		boost::timer::auto_cpu_timer t{"cpu->cpu %ws wall, CPU (%p%)\n"};
//		B_cpu = A_cpu;
//	}
//	{
//		multi::cuda::array<double, 1> A_gpu(size(A_cpu));
//		multi::cuda::array<double, 1> B_gpu(size(A_cpu));
//		boost::timer::auto_cpu_timer t{"gpu->gpu %ws wall, CPU (%p%)\n"};
//		B_gpu = A_gpu;
//		cudaDeviceSynchronize();
//	}
//	{
//		multi::cuda::array<double, 1> B_gpu(size(A_cpu));
//		boost::timer::auto_cpu_timer t{"cpu->gpu %ws wall, CPU (%p%)\n"};
//		B_gpu = A_cpu;
//	}
//	multi::cuda::array<double, 1> C_cpu(size(A_cpu));
//	{
//		multi::cuda::array<double, 1> B_gpu(size(A_cpu));
//		boost::timer::auto_cpu_timer t{"gpu->cpu %ws wall, CPU (%p%)\n"};
//		C_cpu = B_gpu;
//		cudaDeviceSynchronize();
//	}
//	{
//		multi::cuda::managed::array<double, 1> const A_mng(size(A_cpu));
//		multi::cuda::managed::array<double, 1> B_mng(size(A_cpu));
//		{
//			boost::timer::auto_cpu_timer t{"cold mng->mng %ws wall, CPU (%p%)\n"};
//			B_mng = A_mng;
//			cudaDeviceSynchronize();
//		}
//		{
//			boost::timer::auto_cpu_timer t{"haut mng->mng %ws wall, CPU (%p%)\n"};
//			B_mng = A_mng;
//			cudaDeviceSynchronize();
//		}
//	}

////	multi::array<double, 1> Bcpu(3); 
////	Bcpu = B;
////	BOOST_REQUIRE( Bcpu[1] == 3. );
//}

//#if 0
//BOOST_AUTO_TEST_CASE(multi_memory_adaptors_cuda_copy_2D){
//	multi::array<double, 1> A(50, 99.);
//	multi::cuda::array<double, 1> B(50);
//	BOOST_REQUIRE( size(B) == 50 );

////	using std::copy_n;
////	using std::copy;
//	using boost::multi::adl::copy_n;
////	copy_n(&A[0], size(A), &B[0]);
//	copy_n(begin(A), size(A), begin(B));

//	multi::cuda::array<double, 1> D(50);
//	copy_n(begin(B), size(B), begin(D));

//	multi::array<double, 1> C(50, 88.);
//	copy_n(begin(D), size(D), begin(C));
////	C = B;

////	BOOST_REQUIRE( C == A );
//}

//BOOST_AUTO_TEST_CASE(multi_cuda_managed_array_initialization_complex){
//	multi::cuda::managed::array<complex, 1> B = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
//	multi::array<complex, 1> Bcpu(3); 
//	Bcpu = B;
//	BOOST_REQUIRE( Bcpu[1] == 3. + 1.*I );
//}

//namespace utf = boost::unit_test;

//#if 0
//BOOST_AUTO_TEST_CASE(multi_memory_adaptors_cuda_algorithm, *utf::disabled()){
//	BOOST_REQUIRE(false);
//	multi::cuda::array<double, 1> const A(10, 99.);
//	multi::cuda::array<double, 1> B(10, 88.);
//	B = A;

//	B() = A();

////	B.assign({1., 2., 3., 4.});
//	B = {1., 2., 3., 4.};
//	BOOST_REQUIRE( size(B) == 4 );

//	B().assign({11., 22., 33., 44.});//.begin(), il.end());
////	BOOST_REQUIRE( B[2] == 33. );
/////	B.assign


////	multi::cuda::array<double, 2> B({10, 10}, 88.);
////	B = A;
//#if 0
//	{
//		cuda::allocator<double> calloc;
//		std::size_t n = 2e9/sizeof(double);
//		[[maybe_unused]] cuda::ptr<double> p = calloc.allocate(n);
//		{
//			boost::timer::auto_cpu_timer t;
//		//	using std::fill_n; fill_n(p, n, 99.);
//		}
//	//	assert( p[0] == 99. );
//	//	assert( p[n/2] == 99. );
//	//	assert( p[n-1] == 99. );
//		[[maybe_unused]] cuda::ptr<double> q = calloc.allocate(n);
//		{
//			boost::timer::auto_cpu_timer t;
//		//	multi::omp_copy_n(static_cast<double*>(p), n, static_cast<double*>(q));
//		//	using std::copy_n; copy_n(p, n, q);
//			using std::copy; copy(p, p + n, q);
//		}
//		{
//			boost::timer::auto_cpu_timer t;
//		//	using std::copy_n; copy_n(p, n, q);
//		}
//	//	assert( q[23] == 99. );
//	//	assert( q[99] == 99. );
//	//	assert( q[n-1] == 99. );
//	}
//	{
//		multi::array<double, 1> const A(100);//, double{99.});
//		multi::array<double, 1, cuda::allocator<double>> A_gpu = A;
//				#pragma GCC diagnostic push // allow cuda element access
//				#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
//	//	assert( A_gpu[1] == A_gpu[0] );
//				#pragma GCC diagnostic pop
//	}
//#endif
//	{
//		multi::array<double, 2> A({32, 8}, 99.);
//		multi::array<double, 2, cuda::allocator<double>> A_gpu({32, 8}, 0.);// = A;//({32, 8000});// = A;
//	}
//	
//}
//#endif
//#endif
//#endif
#endif
