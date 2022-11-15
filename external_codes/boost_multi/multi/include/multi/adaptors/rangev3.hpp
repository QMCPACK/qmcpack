#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&c++ -D_TEST_MULTI_ADAPTORS_RANGEV3 -Wall -Wfatal-errors $0.cpp -o$0x &&$0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2019

// This header is to make some exotic cases of Multi iterators with Ranges v3
#include "../array.hpp"

#include<range/v3/range_fwd.hpp>

#ifndef MULTI_ADAPTORS_RANGEV3_HPP
#define MULTI_ADAPTORS_RANGEV3_HPP

namespace ranges{namespace v3{
namespace concepts{ // needed for later version of rangesv3
// this allows to recognize const_iterator as RandomAccessIterator
#if 0
	template<class MA>
	struct common_reference<
		boost::multi::basic_array<typename MA::element, MA::dimensionality, typename MA::element_const_ptr, typename MA::layout_t>&&, 
		MA&
	>{
		using type = boost::multi::basic_array<typename MA::element, MA::dimensionality, typename MA::element_const_ptr, typename MA::layout_t>&&;
	};
	template<class MA>
	struct common_reference<
		MA&, 
		boost::multi::basic_array<typename MA::element, MA::dimensionality, typename MA::element_const_ptr, typename MA::layout_t>&&
	>{
		using type = boost::multi::basic_array<typename MA::element, MA::dimensionality, typename MA::element_const_ptr, typename MA::layout_t>&&;
	};
#endif
}
}}

#ifdef _TEST_MULTI_ADAPTORS_RANGEV3

#include <range/v3/all.hpp>

namespace multi = boost::multi;

int main(){

	using I = multi::array<double, 3>::const_iterator;
	static_assert( ranges::RandomAccessIterator<I>{}, "!");

}
#endif

#endif

