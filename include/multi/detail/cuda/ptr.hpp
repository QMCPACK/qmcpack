#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0.cpp) && nvcc --compiler-options -std=c++14,-Wall,-Wextra,-Wpedantic`#,-Wfatal-errors` -D_TEST_BOOST_MULTI_DETAIL_MEMORY_CUDA_PTR $0.cpp -o $0x && $0x && rm $0x $0.cpp; exit
#endif

#ifndef BOOST_MULTI_DETAIL_MEMORY_CUDA_PTR_HPP
#define BOOST_MULTI_DETAIL_MEMORY_CUDA_PTR_HPP

#include<cuda_runtime.h> // cudaError_t

#include<cassert>
#include<cstddef> // nullptr_t
#include<iterator> // random_access_iterator_tag

#include<type_traits> // is_const

namespace boost {
namespace multi {namespace detail {
namespace memory {namespace cuda {

template<class T> struct ref;

//template<class T> class allocator;

template<typename T, typename Ptr = T*> struct ptr;

template<typename T, typename Ptr>
class ptr {
 protected:
	using impl_t = Ptr;
	impl_t impl_;

 private:
	template<class TT> friend class allocator;
	template<typename, typename> friend class ptr;
	template<class TT, typename = std::enable_if_t<not std::is_const<TT>{}>> 
	ptr(ptr<TT const> const& p) : impl_{const_cast<T*>(impl_)}{}
	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);

 public:
	explicit ptr(impl_t impl) : impl_{impl}{}
	ptr() = default;
	ptr(ptr const&) = default;

	// cppcheck-suppress noExplicitConstructor ; initialize from nullptr
	ptr(std::nullptr_t n) : impl_{n} {}

	template<class Other, typename = decltype(impl_t{std::declval<Other const&>().impl_})>
	// cppcheck-suppress noExplicitConstructor ; TODO(correaa) : implement implicit propagation
	ptr(Other const& o) : impl_{o.impl_}{}

	ptr& operator=(ptr const&) = default;
	auto operator==(ptr const& other) const{return impl_==other.impl_;}
	auto operator!=(ptr const& other) const{return impl_!=other.impl_;}

	using element_type = typename std::pointer_traits<impl_t>::element_type;
	using difference_type = typename std::pointer_traits<impl_t>::difference_type;
	using value_type = T;
	using pointer = ptr<T>;
	using reference = ref<element_type>;
	using iterator_category = typename std::iterator_traits<impl_t>::iterator_category;
//	using iterator_concept  = typename std::iterator_traits<impl_t>::iterator_concept;
	explicit operator bool() const{return impl_;}

	ptr& operator++(){++impl_; return *this;}
	ptr& operator--(){--impl_; return *this;}
	ptr  operator++(int){auto tmp = *this; ++(*this); return tmp;}
	ptr  operator--(int){auto tmp = *this; --(*this); return tmp;}
	ptr& operator+=(typename ptr::difference_type n){impl_+=n; return *this;}
	ptr& operator-=(typename ptr::difference_type n){impl_+=n; return *this;}
	ptr operator+(typename ptr::difference_type n) const{return ptr{impl_ + n};}
	ptr operator-(typename ptr::difference_type n) const{return ptr{impl_ - n};}
	ref<element_type> operator*() const{return {*this};}
	ref<element_type> operator[](difference_type n){return *operator+(n);}//*((*this)+n);} 
	friend ptr to_address(ptr const& p){return p;}
	typename ptr::difference_type operator-(ptr const& other) const{return impl_-other.impl_;}
};

template<typename Ptr>
class ptr<void const, Ptr>{
	using T = void const;
	using impl_t = Ptr;
	impl_t impl_;
	template<typename, typename> friend class ptr;
	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);
	explicit ptr(impl_t impl) : impl_{impl}{}

 public:
	ptr() = default;
	ptr(ptr const&) = default;

	// cppcheck-suppress noExplicitConstructor ; initialize from nullptr
	ptr(std::nullptr_t n) : impl_{n}{}

	template<class Other, typename = decltype(impl_t{std::declval<Other const&>().impl_})>
	// cppcheck-suppress noExplicitConstructor ; any pointer is convertible to void pointer
	ptr(Other const& o) : impl_{o.impl_}{}

	ptr& operator=(ptr const&) = default;

	using pointer = ptr<T>;
	using element_type = typename std::pointer_traits<impl_t>::element_type;
	using difference_type = void;//typename std::pointer_traits<impl_t>::difference_type;
	explicit operator bool() const{return impl_;}
	auto operator==(ptr const& other) const{return impl_==other.impl_;}
	auto operator!=(ptr const& other) const{return impl_!=other.impl_;}
	friend ptr to_address(ptr const& p){return p;}
};

template<typename Ptr>
struct ptr<void, Ptr>{
 protected:
	using T = void;
	using impl_t = Ptr;
	impl_t impl_;

 private:
	explicit ptr(impl_t impl) : impl_{impl} {}
	ptr(ptr<void const> const& p) : impl_{const_cast<void*>(p.impl_)} {}
	template<class TT> friend ptr<TT> const_pointer_cast(ptr<TT const> const&);
	template<class, class> friend class ptr;

 public:
	ptr() = default;
	ptr(ptr const&) = default;

	// cppcheck-suppress noExplicitConstructor ; initialize from nullptr
	ptr(std::nullptr_t n) : impl_{n} {}

	template<class Other, typename = decltype(impl_t{std::declval<Other const&>().impl_})>
	// cppcheck-suppress noExplicitConstructor ; any pointer is convertible to void pointer
	ptr(Other const& o) : impl_{o.impl_} {}

	ptr& operator=(ptr const&) = default;

	auto operator==(ptr const& other) const{return impl_==other.impl_;}
	auto operator!=(ptr const& other) const{return impl_!=other.impl_;}

	using pointer = ptr<T>;
	using element_type = typename std::pointer_traits<impl_t>::element_type;
	using difference_type = void;// typename std::pointer_traits<impl_t>::difference_type;

	explicit operator bool() const{return impl_;}
	friend ptr to_address(ptr const& p){return p;}
};

template<class T> 
ptr<T> const_pointer_cast(ptr<T const> const& p){return {p.impl_};}

template<class... Fs> struct overload{}; //template<> struct overload<>{};
template<class F, class... Fs> 
struct overload<F, Fs...> : F, Fs...{
	overload(F f, Fs... fs) : F{std::move(f)}, Fs{std::move(fs)}...{}
	using F::operator();
};
template<class... Fs> 
overload<Fs...> make_overload(Fs&&... fs){return {std::forward<Fs>(fs)...};}

template<class T> struct ref;

template<> struct ref<void>{};

template<class T>
struct ref : private ptr<T>{
	using value_type = T;
	using reference = value_type&;
	using pointer = ptr<T>;

private:
	explicit ref(pointer p) : ptr<T>{std::move(p)}{}
	friend class ptr<T>;
	ptr<T> operator&(){return *this;}
	struct skeleton_t {
		std::array<char, sizeof(T)> buff;  // char buff[sizeof(T)];
		T* p_;
		explicit skeleton_t(T* p) : p_{p} {cudaError_t s = cudaMemcpy(buff.data(), p_, buff.size(), cudaMemcpyDeviceToHost); assert(s == cudaSuccess);}
		operator T&() && {return reinterpret_cast<T&>(buff);}
		void conditional_copyback_if_not(std::false_type) const {
			cudaError_t s = cudaMemcpy(p_, buff.data(), buff.size(), cudaMemcpyHostToDevice); assert(s == cudaSuccess);
		}
		void conditional_copyback_if_not(std::true_type) const {}
		~skeleton_t(){conditional_copyback_if_not(std::is_const<T>{});}
	};
	skeleton_t skeleton()&&{return {this->impl_};}

 public:
	// cppcheck-suppress noExplicitConstructor ; bug in cppcheck 2.3
	ref(ref&& r) : ptr<T>{r}{}
	ref& operator=(ref const&)& = delete;

 private:
	ref& move_assign(ref&& other, std::true_type)&{
		cudaError_t s = cudaMemcpy(this->impl_, other.impl_, sizeof(T), cudaMemcpyDeviceToDevice); assert(s == cudaSuccess);
		return *this;
	}
	ref& move_assign(ref&& other, std::false_type)&{
		cudaError_t s = cudaMemcpy(this->impl_, other.impl_, sizeof(T), cudaMemcpyDeviceToDevice); assert(s == cudaSuccess);
		return *this;
	}
public:
	ref&& operator=(ref&& other)&&{return std::move(move_assign(std::move(other), std::is_trivially_copy_assignable<T>{}));}
private:
public:
	template<class Other>
	auto operator+(Other&& o)&&
	->decltype(std::move(*this).skeleton() + std::forward<Other>(o)) {
		return std::move(*this).skeleton() + std::forward<Other>(o); }
//	template<class Self, class O, typename = std::enable_if_t<std::is_same<std::decay_t<Self>, ref>{}> > 
//	friend auto operator+(Self&& self, O&& o)
//	->decltype(std::forward<Self>(self).skeleton() + std::forward<O>(o)){
//		return std::forward<Self>(self).skeleton() + std::forward<O>(o);}
	ref&& operator=(value_type const& t) && {
		make_overload(
			[&](std::true_type ) {cudaError_t s= cudaMemcpy(this->impl_, std::addressof(t), sizeof(T), cudaMemcpyHostToDevice); assert(s == cudaSuccess);},
			[&](std::false_type) {
				std::array<char, sizeof(T)> buff;  // char buff[sizeof(T)];
				cudaError_t s1 = cudaMemcpy(buff.data(), this->impl_, buff.size(), cudaMemcpyDeviceToHost); assert(s1 == cudaSuccess);
				reinterpret_cast<T&>(buff) = t;
				cudaError_t s2 = cudaMemcpy(this->impl_, buff.data(), buff.size(), cudaMemcpyHostToDevice); assert(s2 == cudaSuccess);
			}
		)(std::is_trivially_copy_assignable<T>{});
		return std::move(*this);
	}
	template<class Other>
	decltype(auto) operator==(ref<Other>&& other) && {
		std::array<char, sizeof(T)> buff1;  // char buff1[sizeof(T)];
		cudaError_t s1 = cudaMemcpy(buff1.data(), this->impl_, buff1.size(), cudaMemcpyDeviceToHost); assert(s1 == cudaSuccess);
		std::array<char, sizeof(Other)> buff2;  // char buff2[sizeof(Other)];
		cudaError_t s2 = cudaMemcpy(buff2.data(), other.impl_, buff2.size(), cudaMemcpyDeviceToHost); assert(s2 == cudaSuccess);
		return reinterpret_cast<T const&>(buff1) == reinterpret_cast<Other const&>(buff2);
	}
	template<class Other>
	decltype(auto) operator!=(ref<Other>&& other) && {
		std::array<char, sizeof(T)> buff1;  // char buff1[sizeof(T)];
		cudaError_t s1 = cudaMemcpy(buff1.data(), this->impl_, buff1.size(), cudaMemcpyDeviceToHost); assert(s1 == cudaSuccess);
		std::array<char, sizeof(Other)> buff2;  // char buff2[sizeof(Other)];
		cudaError_t s2 = cudaMemcpy(buff2.data(), other.impl_, buff2.size(), cudaMemcpyDeviceToHost); assert(s2 == cudaSuccess);
		return reinterpret_cast<T const&>(buff1) != reinterpret_cast<Other const&>(buff2);
	}
	operator T() && {
		static_assert(not std::is_same<T, void>{}, "!");
		std::array<char, sizeof(T)> buff;  // char buff[sizeof(T)];
		cudaError_t s = cudaMemcpy(buff.data(), this->impl_, buff.size(), cudaMemcpyDeviceToHost); assert(s == cudaSuccess );
		return std::move(reinterpret_cast<T&>(buff));
	}
	template<class Other, typename = decltype(std::declval<T&>()+=std::declval<Other&&>())>
	decltype(auto) operator+=(Other&& o) && {std::move(*this).skeleton()+=o; return *this;}
	template<class Other, typename = decltype(std::declval<T&>()-=std::declval<Other&&>())>
	decltype(auto) operator-=(Other&& o) && {std::move(*this).skeleton()-=o;}
	friend void swap(ref&& a, ref&& b) {T tmp = std::move(a); a = std::move(b); b = std::move(tmp);}
	decltype(auto) operator++() && {++(std::move(*this).skeleton()); return *this;}
	decltype(auto) operator--() && {--(std::move(*this).skeleton()); return *this;}
};

}  // end namespace cuda
}  // end namespace memory
}  // end namespace detail

}  // end namespace multi
}  // end namespace boost

#ifdef _TEST_BOOST_MULTI_DETAIL_MEMORY_CUDA_PTR

#include<memory>
#include<iostream>
#include "../../../multi/array.hpp"
#include "../cuda/allocator.hpp"

namespace boost{
namespace multi{namespace cuda{
	template<class T, dimensionality_type D>
	using array = multi::array<T, D, multi::detail::memory::cuda::allocator<T>>;
}}
}

namespace multi = boost::multi;
namespace cuda = multi::detail::memory::cuda;

void add_one(double& d){d += 1.;}
template<class T>
void add_one(T&& t){std::forward<T>(t) += 1.;}

int main() {

	static_assert(std::is_same<typename std::pointer_traits<cuda::ptr<double>>::element_type, double>{}, "!");
	cuda::allocator<double> calloc;
	cuda::ptr<double> p = calloc.allocate(100);
	cuda::ptr<void> v = p;
	cuda::ptr<void const> vc{v};
	v = const_pointer_cast(vc);
	assert( vc == v );
	std::pointer_traits<decltype(p)>::rebind<double const> pc = p; // cuda::ptr<double const> pc = p;
	assert( pc == p );
	using cuda::const_pointer_cast;
	auto end = p + 100;
	auto rbegin = std::make_reverse_iterator(end);
	auto rend = std::make_reverse_iterator(p);
	std::transform(rbegin, rend, rbegin, [](auto&& e){return std::forward<decltype(e)>(e) + 99.;});
	assert( p[11] == 99. );
	p[33] = 123.;
	p[99] = 321.;
//	p[33] += 1;
	add_one(p[33]);
	double p33 = p[33];
	assert( p33 == 124. );
	assert( p[33] == 124. );
	assert( p[33] == p[33] );
	swap(p[33], p[99]);
	assert( p[99] == 124. );
	assert( p[33] == 321. );
	std::cout << p[33] << std::endl;
	calloc.deallocate(p, 100);

	multi::array<double, 1, cuda::allocator<double>> arr2(multi::array<double, 1>::extensions_type{100l}, 999.);

	assert(size(arr2) == 100);
}
#endif
#endif
