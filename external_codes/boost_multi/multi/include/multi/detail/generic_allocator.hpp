// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2021 Alfredo A. Correa

#ifndef BOOST_MULTI_DETAIL_GENERIC_ALLOCATOR_HPP
#define BOOST_MULTI_DETAIL_GENERIC_ALLOCATOR_HPP

#include "../detail/memory.hpp"

#include<cassert>
#include<memory>

#if __cplusplus > 201703L
#include<memory_resource>
#endif
// static_assert(__cpp_lib_experimental_memory_resources==201402, "!");
#include <utility>  // for forward

namespace boost {
namespace multi {

template<class MR>
auto allocator_of(MR& mr)
->decltype(mr->allocator()) {
	return mr->allocator(); }

inline std::allocator<char>& allocator_of(...) {
	static std::allocator<char> instance;
	return instance;
}

template<class T, class MemoryResource
#if(__cpp_lib_memory_resource==201603L)
	= std::pmr::memory_resource
#endif
>
class generic_allocator{
	using memory_resource_type = MemoryResource;
	memory_resource_type* mr_;
	template<class TT, class MR2> friend class generic_allocator;

 public:
	using value_type = T;
	using pointer = typename std::pointer_traits<decltype(std::declval<MemoryResource*>()->allocate(0))>::template rebind<value_type>;
	using difference_type = typename std::pointer_traits<pointer>::difference_type;
	using size_type =  std::make_unsigned_t<difference_type>;

	generic_allocator() : mr_{nullptr} {}

	// cppcheck-suppress noExplicitConstructor ; allocators are pointers to memory resources
	generic_allocator(memory_resource_type* mr) : mr_{mr} {}  // NOLINT(runtime/explicit)

	template<typename T2>
	generic_allocator(generic_allocator<T2, MemoryResource> const& other)
	: mr_{other.mr_} {}

	bool operator==(generic_allocator const& o) const {return mr_ == o.mr_;}
	bool operator!=(generic_allocator const& o) const {return not(o==*this);}

	pointer allocate(size_type n) {
		if(n and !mr_) throw std::bad_alloc{};
		return static_cast<pointer>(mr_->allocate(n*sizeof(value_type)));
	}
	void deallocate(pointer p, size_type n) {
		if(n==0 and p == nullptr) return;
		mr_->deallocate(p, n*sizeof(value_type));
	}
	template<class... Args>
	void construct(pointer p, Args&&... args) {
//  ->decltype(allocator_traits<std::decay_t<decltype(allocator_of(std::declval<memory_resource_type&>()))>>::construct(allocator_of(*mr_), p, std::forward<Args>(args)...)){
	//  mr_->allocator().construct(p, std::forward<Args>(args)...);
	//  using TA = allocator_traits<std::decay_t<decltype(allocator_of(mr_))>>;
		allocator_traits<std::decay_t<decltype(allocator_of(mr_))>>::construct(allocator_of(mr_), p, std::forward<Args>(args)...);
	}
	decltype(auto) destroy(pointer p){
	//  mr_->allocator().destroy(p);
		allocator_traits<std::decay_t<decltype(allocator_of(mr_))>>::destroy(allocator_of(mr_), p);
	}
};

}  // end namespace multi
}  // end namespace boost

#ifdef _TEST_BOOST_MULTI_DETAIL_GENERIC_ALLOCATOR

#include<vector>
#include "../array.hpp"
#include<iostream>

namespace multi = boost::multi;
using std::cout;

int main() {
#if 1
	multi::generic_allocator<double, std::pmr::memory_resource> ga(std::pmr::get_default_resource());
	double* p = ga.allocate(1);
	std::allocator_traits<multi::generic_allocator<double, std::pmr::memory_resource>>::construct(ga, p, 8.);
//  ga.construct(p, 8.);
	assert( *p == 8. );

	std::vector<double, multi::generic_allocator<double, std::pmr::memory_resource>> v(100, std::pmr::get_default_resource());
//  std::vector v(100, 1.2, multi::allocator<double>{}); // needs C++17 CTAD
	multi::array<double, 2, multi::generic_allocator<double, std::pmr::memory_resource>> m({2, 4}, 0., std::pmr::get_default_resource());
//  multi::array m({2,4}, 0., pmr::get_default_resource()); // needs C++17 CTAD
	m[1][3] = 99.;
	assert( m[1][3] == 99. );
#endif
}
#endif
#endif

