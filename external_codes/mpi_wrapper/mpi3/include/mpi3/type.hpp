// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#ifndef BOOST_MPI3_TYPE_HPP
#define BOOST_MPI3_TYPE_HPP

#include <mpi3/core.hpp>
#include <mpi3/detail/datatype.hpp>

#if defined(__NVCC__)
#include <thrust/complex.h>
#endif

#include <complex>
#include <functional>  // for std::invoke
#include <map>
#include <typeindex>

namespace boost {
namespace mpi3 {

template<
	class MultiIt,
	class Size                                            = typename MultiIt::difference_type,
	class Stride                                          = typename MultiIt::stride_type,
	std::enable_if_t<(MultiIt::dimensionality >= 1), int> = 0,
	typename Element = typename MultiIt::element, typename DataType = detail::basic_datatype<Element>,
	std::enable_if_t<detail::is_basic<Element>{}, int> = 0>
typename MultiIt::element_ptr base(MultiIt first) { return first.base(); }

template<class T> T* base(T* p) { return p; }

struct committed_type {
 private:
	MPI_Datatype impl_;
	explicit committed_type(MPI_Datatype dt) : impl_{dt} {}
	friend struct type;

 public:
	auto get() const -> MPI_Datatype {return impl_;}
	explicit operator MPI_Datatype() const { return impl_; }
};

struct type {
	explicit type(MPI_Datatype const& dt) noexcept : impl_{dt} {  // NOLINT(bugprone-exception-escape) TODO(correaa) improve this global initialization
	 	if(mpi3::initialized()) {  // cppcheck-suppress[throwInNoexceptFunction]; TODO(correaa) improve this global initialization
	 		MPI_(Type_dup)(dt, &impl_);
	 	}
	}

	template<class T>
	explicit type(detail::basic_datatype<T> bd) : impl_(bd) {
	}

	// template<class T, typename = decltype(detail::basic_datatype<T>::value_f())>
	// explicit type(T const* /*unused*/) { MPI_Type_dup(detail::basic_datatype<T>::value_f(), &impl_); }

	template<
		class T,
		std::enable_if_t<not std::is_same<T*, MPI_Datatype>{}, int>                                   = 0,
		std::enable_if_t<std::is_trivially_copy_assignable<T>{} and (not detail::is_basic<T>{}), int> = 0>
	explicit type(T const* /*type*/)
	: type{type{MPI_BYTE}.contiguous(sizeof(T))} {}

	template<
		class MultiIt, class Size = typename MultiIt::difference_type, class Stride = typename MultiIt::stride_type, std::enable_if_t<MultiIt::dimensionality == 1, int> = 0,
		typename E = typename MultiIt::element, typename = decltype(detail::basic_datatype<E>::value_f()),
		std::enable_if_t<detail::is_basic<E>{}, int> = 0>
	explicit type(MultiIt first) : type{type{first.base()}.vector(1, 1, first.stride() * sizeof(E)).resize(0, first.stride() * sizeof(E))} {}

 private:
	template<class F, class Tuple, std::size_t... I>
	static decltype(auto) apply_impl(F&& f, Tuple const& t, std::index_sequence<I...> /*012*/) {
		return std::forward<F>(f)(std::get<I>(t)...);
	}
	template<class F, class Tuple>
	static decltype(auto) apply(F&& f, Tuple const& t) {
		return apply_impl(
			std::forward<F>(f), std::forward<Tuple>(t),
			std::make_index_sequence<std::tuple_size<Tuple>{}>{}
		);
	}

 public:
	template<
		class MultiIt, class Stride = typename MultiIt::stride_type,
		std::size_t D = MultiIt::dimensionality, std::enable_if_t<(D >= 2), int> = 0,
		typename E = typename MultiIt::element, typename = decltype(detail::basic_datatype<E>::value_f()),
		std::enable_if_t<detail::is_basic<E>{}, int> = 0>
	explicit type(MultiIt first) : type{first.base()} {
		auto const strides = apply([](auto... e) { return std::array<Stride, D - 1>{static_cast<Stride>(e)...}; }, first->strides());  // NOLINT(altera-id-dependent-backward-branch) TODO(correaa) investigate
		auto const sizes   = apply([](auto... e) { return std::array<Stride, D - 1>{static_cast<Stride>(e)...}; }, first->sizes());  // NOLINT(altera-id-dependent-backward-branch) TODO(correaa) investigate
		for(Stride i = 1; i != Stride{strides.size()} + 1; ++i) {  // NOLINT(altera-id-dependent-backward-branch,altera-unroll-loops) TODO(correaa) use an algorithm
			(*this) = this->vector(sizes[sizes.size() - i], 1, strides[strides.size() - i]).resize(0, strides[strides.size() - i] * sizeof(E));
		}
		(*this) = this->vector(1, 1, first.stride() * sizeof(E)).resize(0, first.stride() * sizeof(E));
	}

	MPI_Datatype impl_ = MPI_DATATYPE_NULL;  // NOLINT(misc-non-private-member-variables-in-classes) TODO(correaa)

	type() = default;// {std::clog << "ctor type()" << std::endl;}
	type(type const& other) { MPI_Type_dup(other.impl_, &impl_); }
	type(type&& other) noexcept : impl_{std::exchange(other.impl_, MPI_DATATYPE_NULL)} {}  // TODO(correaa) consider not making it default constructible or movable

	type& operator=(type const& other) {
		type tmp(other);
		swap(tmp);
		return *this;
	}
	type& operator=(type&& other) noexcept {  // TODO(correaa) consider not making it default constructible or movable
		type tmp(std::move(other));
		swap(tmp);
		return *this;
	}

	void     swap(type& other) { std::swap(impl_, other.impl_); }
	explicit operator MPI_Datatype() const& {
		MPI_Type_commit(const_cast<MPI_Datatype*>(&impl_));  // NOLINT(cppcoreguidelines-pro-type-const-cast) TODO(correaa)
		return impl_;
	}

	auto get() const -> MPI_Datatype {return impl_;}

	auto operator&() const& -> type const* { return this; }  // NOLINT(google-runtime-operator)
	auto operator&() && -> MPI_Datatype {  // NOLINT(google-runtime-operator)
		MPI_Type_commit(const_cast<MPI_Datatype*>(&impl_));  // NOLINT(cppcoreguidelines-pro-type-const-cast) TODO(correaa)
		return impl_;
	}
	auto operator&() & -> MPI_Datatype {  // NOLINT(google-runtime-operator)
		MPI_Type_commit(const_cast<MPI_Datatype*>(&impl_));  // NOLINT(cppcoreguidelines-pro-type-const-cast) TODO(correaa)
		return impl_;
	}

	committed_type commit() && {
		MPI_Type_commit(const_cast<MPI_Datatype*>(&impl_));  // NOLINT(cppcoreguidelines-pro-type-const-cast) TODO(correaa)
		return committed_type{std::exchange(impl_, MPI_DATATYPE_NULL)};
	}
	template<class T> void commit_as(T const& /*unused*/) { return commit_as<T>(); }
	~type() noexcept {
		// std::clog << "dtor type" << std::endl;
		try {
			if(mpi3::initialized() and not mpi3::finalized()) {  // TODO(correaa) types defined statically will generally leak on exit
				if(impl_ != MPI_DATATYPE_NULL) {
					MPI_Type_free(&impl_);
				}
			}
		} catch(...) {
		}
	}

	type contiguous(int count) const {
		type      ret;  // NOLINT() delayed init
		int const s = MPI_Type_contiguous(count, impl_, &ret.impl_);  // TODO(correaa) modernize calls
		if(s != MPI_SUCCESS) {
			throw std::runtime_error{"cannot build contiguous type"};
		}
		ret.set_name("(" + name() + ")[" + std::to_string(count) + "]");
		return ret;
	}
	type vector(int count, int block_length, int stride) const {  // element units, hvector for bytes
		type ret;
		MPI_Type_vector(count, block_length, stride, impl_, &ret.impl_);
		using std::to_string;
		ret.set_name("(" + name() + ")[" + to_string(count) + "," + to_string(block_length) + ":" + to_string(stride) + "]");
		return ret;
	}
	type resize(MPI_Aint lower_bound, MPI_Aint extent) const {
		type ret;
		MPI_Type_create_resized(impl_, lower_bound, extent, &ret.impl_);
		return ret;
	}
	type stride(int stride) const { return resize(0, static_cast<MPI_Aint>(stride) * size()); }

	// MPI_Type_struct is deprecated
	static type struct_(std::initializer_list<type> il) {  // NOLINT(readability-identifier-naming) meta
		type                  ret;
		std::vector<int>      blocklen(il.size(), 1);
		std::vector<MPI_Aint> disp;
		disp.reserve(il.size());
		//	std::vector<MPI_Datatype> array_of_types;
		//	array_of_types.reserve(il.size());
		MPI_Aint    current_disp = 0;
		std::string new_name     = "{";
		std::for_each(il.begin(), il.end(), [&il, &disp, &current_disp, &new_name](auto const& e) {
			disp.push_back(current_disp);
			current_disp += e.size();
			new_name += (&e != il.begin() ? ", " : "") + e.name();
			// array_of_types.push_back(e.impl_);
		});

		MPI_Type_create_struct(
			static_cast<int>(il.size()),
			blocklen.data(),
			disp.data(),
			&il.begin()->impl_,
			&ret.impl_
		);

		ret.name(new_name);
		return ret;
	}

	type operator[](int count) const { return contiguous(count); }
	type operator()(int stride) const {
		//	assert( stride == 2 );
		return vector(1, 1, stride);
	}
	int size() const {
		int ret;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
		MPI_Type_size(impl_, &ret);
		return ret;
	}

	std::string name() const {
		std::array<char, MPI_MAX_OBJECT_NAME> name{};
		int                                   namelen;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		MPI_Type_get_name(impl_, name.data(), &namelen);
		return {name.data(), static_cast<std::size_t>(namelen)};
	}
	void name(std::string const& s) { set_name(s); }
	void set_name(std::string const& s) { MPI_Type_set_name(impl_, s.c_str()); }  // NOLINT(readability-make-member-function-const) this is not really const

	MPI_Aint extent() const {
		MPI_Aint  lb;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		MPI_Aint  ext;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		int const s = MPI_Type_get_extent(impl_, &lb, &ext);  // TODO(correaa) modernize calls
		if(s != MPI_SUCCESS) {
			throw std::runtime_error{"cannot extent"};
		}
		return ext;
	}
	MPI_Aint lower_bound() const {
		MPI_Aint  lb;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		MPI_Aint  ext;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		int const s = MPI_Type_get_extent(impl_, &lb, &ext);  // TODO(correaa) modernize calls
		if(s != MPI_SUCCESS) {
			throw std::runtime_error{"cannot lower bound"};
		}
		return lb;
	}
	MPI_Aint upper_bound() const {
		MPI_Aint  lb;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		MPI_Aint  ext;  // NOLINT(cppcoreguidelines-init-variables) delayed init
		int const s = MPI_Type_get_extent(impl_, &lb, &ext);  // TODO(correaa) modernize call
		if(s != MPI_SUCCESS) {
			throw std::runtime_error{"cannot lower bound"};
		}
		return lb + ext;
	}
	type operator,(type const& other) const {
		type                        ret;
		int const                   count          = 2;
		std::array<int, 2>          blocklen       = {1, 1};
		std::array<MPI_Aint, 2>     disp           = {0, this->size()};
		std::array<MPI_Datatype, 2> array_of_types = {impl_, other.impl_};
		MPI_Type_create_struct(count, blocklen.data(), disp.data(), array_of_types.data(), &ret.impl_);
		std::string const newname = name() + ", " + other.name();
		MPI_Type_set_name(ret.impl_, newname.c_str());
		return ret;
	}
	//	static std::map<std::type_index, type const&> registered;
};

// vvv TODO(correaa)
// vvv this will work in clang-tidy 14 https://clang.llvm.org/extra/clang-tidy/
// NOLINTBEGIN(fuchsia-statically-constructed-objects)
static type const  char_{MPI_CHAR};  // NOLINT(fuchsia-statically-constructed-objects)
static type const  unsigned_char{MPI_UNSIGNED_CHAR};
static type const& unsigned_char_ = unsigned_char;  // NOLINT(fuchsia-statically-constructed-objects)
static type const  short_{MPI_SHORT};  // NOLINT(fuchsia-statically-constructed-objects)
static type const  unsigned_short{MPI_UNSIGNED_SHORT};
static type const& unsigned_short_ = unsigned_short;  // NOLINT(fuchsia-statically-constructed-objects)
static type const  int_{MPI_INT};  // NOLINT(fuchsia-statically-constructed-objects)
static type const  unsigned_int_{MPI_UNSIGNED};
static type const& unsigned_int = unsigned_int_;
static type const& unsigned_    = unsigned_int_;  // NOLINT(fuchsia-statically-constructed-objects)
static type const  long_{MPI_LONG};  // NOLINT(fuchsia-statically-constructed-objects)
static type const  unsigned_long{MPI_UNSIGNED_LONG};
static type const& unsigned_long_ = unsigned_long;  // NOLINT(fuchsia-statically-constructed-objects)
static type const  float_{MPI_FLOAT};  // NOLINT(fuchsia-statically-constructed-objects)
static type const  double_{MPI_DOUBLE};  // NOLINT(fuchsia-statically-constructed-objects)
static type const  long_double_{MPI_LONG_DOUBLE};
static type const& long_double = long_double_;  // NOLINT(fuchsia-statically-constructed-objects)
static type const  long_long_int{MPI_LONG_DOUBLE_INT};  // NOLINT(fuchsia-statically-constructed-objects)
static type const  float_int{MPI_FLOAT_INT};  // NOLINT(fuchsia-statically-constructed-objects)
static type const  long_int{MPI_LONG_INT};  // NOLINT(fuchsia-statically-constructed-objects)
static type const  double_int{MPI_DOUBLE_INT};  // NOLINT(fuchsia-statically-constructed-objects)
static type const  short_int{MPI_SHORT_INT};  // NOLINT(fuchsia-statically-constructed-objects)
static type const  int_int{MPI_2INT};
static type const& _2int = int_int;  // NOLINT(fuchsia-statically-constructed-objects)
static type const  long_double_int{MPI_LONG_DOUBLE_INT};  // NOLINT(fuchsia-statically-constructed-objects)
// NOLINTEND(fuchsia-statically-constructed-objects)

template<class T> auto make_type() -> type;

template<> inline auto make_type<char>() -> type { return type{MPI_CHAR}; }
template<> inline auto make_type<unsigned char>() -> type { return type{MPI_UNSIGNED_CHAR}; }

template<> inline auto make_type<std::complex<float>>() -> type { return type{MPI_CXX_FLOAT_COMPLEX}; }
template<> inline auto make_type<std::complex<double>>() -> type { return type{MPI_CXX_DOUBLE_COMPLEX}; }

// template<Complex, class = std::enable_if_t<std::is_same<decltype(Complex{}.real())> > > inline auto make_type<Complex>() -> type { return type{MPI_CXX_FLOAT_COMPLEX}; }
// template<Complex, class = std::enable_if_t<std::is_same<> > inline auto make_type<std::complex<double>>() -> type { return type{MPI_CXX_DOUBLE_COMPLEX}; }

#if defined(__NVCC__)
template<> inline auto make_type<thrust::complex<float>>() -> type { return type{MPI_CXX_FLOAT_COMPLEX}; }
template<> inline auto make_type<thrust::complex<double>>() -> type { return type{MPI_CXX_DOUBLE_COMPLEX}; }
#endif

template<> inline type make_type<double>() { return mpi3::double_; }
template<> inline type make_type<int>() { return mpi3::int_; }

// // TODO(correaa) remove this?
// template<class T> class datatype<T[2]> {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
//  public:
// 	template<class R = decltype(boost::mpi3::type{mpi3::datatype<T>{}()}[2].get())>
// 	R operator()() const {
// 		static auto ret = std::invoke([]() {
// 			assert(boost::mpi3::initialized());
// 			return boost::mpi3::type{mpi3::datatype<T>{}()}[2].commit().get();
// 		});
// 		return ret;
// 	}
// };

// template<class T, std::size_t D> class datatype<std::array<T, D>> {
//  public:
// 	template<class R = decltype(boost::mpi3::type{mpi3::datatype<T>{}()}[D].get())>
// 	R operator()() const {
// 		static auto ret = std::invoke([]() {
// 			assert(boost::mpi3::initialized());
// 			return boost::mpi3::type{mpi3::datatype<T>{}()}[D].commit().get();  // cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays
// 		});
// 		return ret;
// 	}
// };

template<class... Ts> struct struct_;

template<class T1, class T2> struct struct_<std::pair<T1, T2>> : struct_<T1, T2> {};

template<class T1> struct struct_<T1> {
	static auto make() {
		struct dummy {
			T1 t1_;
		};
		const int nitems = 1;
		std::array<int, 1> blocklengths = {1};
		std::array<MPI_Datatype, 1> types = {datatype<T1>{}.get()};
		std::array<MPI_Aint, 1> offsets {
			offsetof(dummy, t1_)
		};
		mpi3::type ret;
		MPI_Type_create_struct(nitems, blocklengths.data(), offsets.data(), types.data(), &ret.impl_);
		return ret;
	}
	auto operator()() const {
		static auto ret = make().commit();
		return ret.get();
	}
};

template<class T1, class T2> struct struct_<T1, T2> {
	static auto make() {
		struct dummy {
			T1 t1_;
			T2 t2_;
		};
		const int nitems = 2;
		std::array<int, 2> blocklengths = {1, 1};
		std::array<MPI_Datatype, 2> types = {datatype<T1>{}.get(), datatype<T2>{}.get()};
		std::array<MPI_Aint, 2> offsets {
			offsetof(dummy, t1_),
			offsetof(dummy, t2_)
		};
		mpi3::type ret;
		MPI_Type_create_struct(nitems, blocklengths.data(), offsets.data(), types.data(), &ret.impl_);
		return ret;
	}
	auto operator()() const {
		static auto ret = make().commit();
		return ret.get();
	}
};

template<class T1, class T2, class T3> struct struct_<T1, T2, T3> {
	static auto make() {
		struct dummy {
			T1 t1_;
			T2 t2_;
			T3 t3_;
		};
		const int nitems = 3;
		std::array<int, 3> blocklengths = {1, 1, 1};
		std::array<MPI_Datatype, 3> types = {
			datatype<T1>{}.get(),
			datatype<T2>{}.get(),
			datatype<T3>{}.get()
		};
		std::array<MPI_Aint, 3> offsets {
			offsetof(dummy, t1_),
			offsetof(dummy, t2_),
			offsetof(dummy, t3_)
		};
		mpi3::type ret;
		MPI_Type_create_struct(nitems, blocklengths.data(), offsets.data(), types.data(), &ret.impl_);
		return ret;
	}
	auto operator()() const {
		static auto ret = make().commit();
		return ret.get();
	}
};

}  // end namespace mpi3
}  // end namespace boost

// #if not __INCLUDE_LEVEL__

// #include "../mpi3/main.hpp"
// #include "../mpi3/communicator.hpp"
// #include "../mpi3/process.hpp"

// #include<iostream>

// namespace mpi3 = boost::mpi3;
// using std::cout;

// struct A{
//	double d[6];
//	long l[20];
// };

// struct B{
//	double d[7];
//	long l[9];
// };

// template<class Archive>
// void serialize(Archive& ar, A& a, const unsigned int){
//	ar & boost::serialization::make_array(a.d, 6);
//	ar & boost::serialization::make_array(a.l, 20);
// }

// int mpi3::main(int argc, char* argv[], mpi3::communicator world){

//	assert(world.size() >= 2);
//	{
//		int value = -1;
//		if(world.rank() == 0){
//			value = 5;
//			world[1] << value; //	world.send_value(value, 1);
//		}else{
//			assert(value == -1);
//			world[0] >> value; //world.receive(&value, &value + 1, 0);
//			assert(value == 5);
//		}
//	}
//	{
//		int buffer[100]; std::fill(&buffer[0], &buffer[100], 0);
//		if(world.rank() == 0){
//			std::iota(&buffer[0], &buffer[100], 0);
//			world[1] << buffer; // world.send_value(buffer, 1);
//		}else{
//			assert(buffer[11]==0);
//			world[0] >> buffer; //	world.receive_value(buffer, 0);
//		//	assert(buffer[11]==11);
//		}
//	}
//	{
//	//	auto Atype = (
//	//		mpi3::double_[6],
//	//		mpi3::long_[20]
//	//	);
//	//	Atype.commit_as<A>();

//		A particle;
//		particle.d[2] = 0.;
//		if(world.rank()==0){
//			particle.d[2] = 5.1;
//			world[1] << particle;
//		}else{
//			assert(particle.d[2]==0.);
//			world[0] >> particle;
//		}
//	}
///*	{
//		auto Btype = (
//			mpi3::double_[6],
//			mpi3::long_[20]
//		);
//		Btype.commit_as<B>();
//		B b;
//		b.d[2] = 0.;
//		if(world.rank()==0){
//			b.d[2] = 5.1;
//			world[1] << b;
//		}else{
//			assert(b.d[2]==0.);
//			world[0] >> b;
//		}
//	}*/
//	return 0;
// #if 0
//	{
//		std::vector<double> v(100);
//		mpi3::type d100 = mpi3::type::double_[100];
//		d100.commit();
//		if(world.rank()==0){
//			v[5] = 6.;
//		//	world.send_n(v.data(), 1, d100);
//		}else{
//		//	world.receive_n(v.data(), 1, d100);
//			assert(v[5] == 6.);
//		}
//	}
// #endif
//	return 0;
//}

// #endif
#endif
