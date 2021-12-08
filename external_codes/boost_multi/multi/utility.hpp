// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#ifndef MULTI_UTILITY_HPP
#define MULTI_UTILITY_HPP

#include "detail/layout.hpp"

#include<memory>  // for allocator<>

#if(__cplusplus >= 201703L)
#include<iterator>  // for std::size (in c++17)
#endif

#include<tuple>  // for tuple<>

namespace boost {
namespace multi {

template<class Array, typename Reference = void, typename Element = void>
struct array_traits;

template<class Array, typename Reference, typename Element>
struct array_traits{
	using reference = typename Array::reference;
	using element   = typename Array::element;
	using element_ptr = typename Array::element_ptr;
	using decay_type = typename Array::decay_type;
	using default_allocator_type = typename Array::default_allocator_type;
};

template<class T, typename = typename T::rank>
       auto has_rank_aux(T const&) -> std::true_type;
inline auto has_rank_aux(...     ) -> std::false_type;

template<class T> struct has_rank : decltype(has_rank_aux(std::declval<T>())){};

template<typename T> struct rank;

template<typename T, typename = std::enable_if_t<has_rank<T>{}> >
constexpr auto rank_aux(T const&) -> typename T::rank;

template<typename T, typename = std::enable_if_t<not has_rank<T>{}> >
constexpr auto rank_aux(T const&) -> std::integral_constant<size_t, std::rank<T>{}>;

template<typename T> struct rank : decltype(rank_aux(std::declval<T>())){};

#if not defined(__cpp_lib_nonmember_container_access) or __cpp_lib_nonmember_container_access < 201411
template<class Container>
constexpr auto size(Container const& con)
-> std::make_signed_t<decltype(con.size())>{
	return static_cast<std::make_signed_t<decltype(con.size())>>(con.size());}
#else
#endif

template<class Pointer, std::enable_if_t<std::is_pointer<Pointer>{}, int> = 0>  // special sfinae trick
constexpr auto stride(Pointer /*unused*/) -> std::ptrdiff_t {return 1;}

template<class Pointer, std::enable_if_t<std::is_pointer<Pointer>{}, int> =0>  // special sfinae trick
constexpr auto base(Pointer d) -> Pointer {return d;}

template<class T, class U>
auto reinterpret_pointer_cast(U* other)
-> decltype(reinterpret_cast<T*>(other)) {return reinterpret_cast<T*>(other);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) : unavoidalbe implementation?

template <class T, std::size_t N>
constexpr auto size(const T(&/*t*/)[N]) noexcept{return multi::size_type{N};}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T, typename = typename T::get_allocator>
       auto has_get_allocator_aux(T const&) -> std::true_type;
inline auto has_get_allocator_aux(...     ) -> std::false_type;

template<class T, std::size_t N>
constexpr auto get_allocator(T(&/*t*/)[N]) noexcept -> std::allocator<std::decay_t<typename std::remove_all_extents<T[N]>::type>>{return {};}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T>
constexpr auto get_allocator(T* const& /*t*/)
->decltype(std::allocator<typename std::iterator_traits<T*>::value_type>{}) {
	return std::allocator<typename std::iterator_traits<T*>::value_type>{}; }

template<class T>
constexpr auto default_allocator_of(T* /*unused*/) {
	return std::allocator<typename std::iterator_traits<T*>::value_type>{};
}

template<class T>
constexpr
auto to_address(T* const& t)
->decltype(t) {
	return t; }

template<class Archive> struct archive_traits {  // TODO(correaa) implement a poors man nvp that works with boost serialization, is it possible?
	template<class T>
	static constexpr auto make_nvp(char const* /*name*/, T& v) -> T&{return v;}
};

template<class T>
auto has_get_allocator_aux(T const& t)->decltype(t.get_allocator(), std::true_type {});

template<class T> struct has_get_allocator : decltype(has_get_allocator_aux(std::declval<T>())){};

template<class T1, class T2, typename Ret = T1>  // std::common_type_t<T1, T2>>
auto common(T1 const& t1, T2 const& t2) -> Ret{
	return t1==t2?t1:Ret{};}

template<class T>
       auto has_num_elements_aux(T const& t)->decltype(std::declval<T const&>().num_elements() + 1, std::true_type {});
inline auto has_num_elements_aux(...       )->decltype(                                             std::false_type{});
template<class T> struct has_num_elements : decltype(has_num_elements_aux(std::declval<T>())){};

template<class A, typename = std::enable_if_t<has_num_elements<A>{}> >
constexpr auto num_elements(A const& arr)
->std::make_signed_t<decltype(arr.num_elements())> {
	return static_cast<std::make_signed_t<decltype(arr.num_elements())>>(arr.num_elements());
}

template<class T>
auto has_size_aux(T const& t)->decltype(t.size(), std::true_type {});
inline auto has_size_aux(...       )->decltype(                  std::false_type{});
template<class T> struct has_size : decltype(has_size_aux(std::declval<T>())){};

template<class T>
auto has_data_elements_aux(T&& t)->decltype(t.data_elements() + 1, std::true_type {});
auto has_data_elements_aux(...  )->decltype(                       std::false_type{});
template<class T> struct has_data_elements : decltype(has_data_elements_aux(std::declval<T>())) {};

template<class T>
auto has_data_aux(T&& t)->decltype(t.data_elements() + 1, std::true_type {});
auto has_data_aux(...  )->decltype(                       std::false_type{});
template<class T> struct has_data : decltype(has_data_aux(std::declval<T>())) {};

template<class A, typename = std::enable_if_t<not has_num_elements<A>{} and has_size<A>{} and has_data<A>{}> >
constexpr auto num_elements(A const& arr)
->std::make_signed_t<decltype(arr.size())> {
	return arr.size();}

template<class A, typename = std::enable_if_t<has_data_elements<A>{}> >
constexpr auto data_elements(A const& arr)
->decltype(arr.data_elements()) {
	return arr.data_elements(); }

template<class T, typename = std::enable_if_t<not std::is_array<T>{}> >
[[deprecated("use constexpr data_elements() or base() to extract pointer")]]
auto data(T& t) {return &t;}

template<class T, typename = std::enable_if_t<not std::is_array<T>{} and not has_data_elements<T>{} and not has_data<T>{}> >
constexpr auto data_elements(T& t) {return &t;}

template<class A> struct num_elements_t: std::integral_constant<std::ptrdiff_t, 1> {};

template<class T, std::size_t N> struct num_elements_t<T[N]>: std::integral_constant<std::ptrdiff_t, (N*num_elements_t<T>{})> {};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T, std::size_t N> struct num_elements_t<T(&)[N]>: num_elements_t<T[N]> {};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T, std::size_t N>
constexpr auto num_elements(const T(&/*t*/)[N]) noexcept {return num_elements_t<T[N]>{};}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class Vector>
constexpr auto num_elements(Vector const& v, std::enable_if_t<std::is_same<typename Vector::pointer, decltype(std::declval<Vector>().data())>{}, int> = 0)  // NOLINT(fuchsia-default-arguments-declarations,hicpp-named-parameter,readability-named-parameter) : special sfinae trick
->std::make_signed_t<decltype(v.size())> {
	return static_cast<std::make_signed_t<decltype(v.size())>>(v.size());}

template<class Vector, typename = std::enable_if_t<std::is_same<typename Vector::pointer, decltype(std::declval<Vector>().data())>{}> >
auto data_elements(Vector const& v)
->decltype(v.data()) {
	return v.data(); }

template <class T, std::size_t N>
constexpr auto stride(const T(&/*unused*/)[N]) noexcept -> std::ptrdiff_t {return num_elements_t<T>{};}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template <class T, std::size_t N>
constexpr auto is_compact(const T(&/*t*/)[N]) noexcept -> bool {return true;}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T, std::size_t N>
constexpr auto offset(const T(&/*t*/)[N]) noexcept -> std::ptrdiff_t {return 0;}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T, std::size_t N>
[[deprecated("use data_elements instead")]]  // this name is bad because when the element belongs to std:: then std::data is picked up by ADL and the
constexpr auto data(T(&t)[N]) noexcept {return data(t[0]);}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T, std::size_t N>
constexpr auto data_elements(T(&t)[N]) noexcept {return data_elements(t[0]);}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T>
       auto has_dimensionality_aux(T const& t)->decltype(T::rank_v, std::true_type {});
inline auto has_dimensionality_aux(...       )->decltype(           std::false_type{});
template<class T> struct has_dimensionality : decltype(has_dimensionality_aux(std::declval<T>())) {};

template<class Container, std::enable_if_t<has_dimensionality<Container>{}, int> =0>
constexpr auto dimensionality(Container const& /*container*/)
->std::decay_t<decltype(Container::rank_v)> {
	return Container::rank_v;}

template<class T>
       auto has_dimensionaliy_member_aux(T const& t)->decltype((size_t(T::rank_v), std::true_type{}));
inline auto has_dimensionaliy_member_aux(...       )->decltype(                    std::false_type{});
template<class T> struct has_dimensionality_member : decltype(has_dimensionaliy_member_aux(std::declval<T>())){};

template<class T, typename = std::enable_if_t<not has_dimensionality_member<T>{}>>
constexpr auto dimensionality(T const&/*, void* = nullptr*/) {return 0;}

template<class T, std::size_t N>
constexpr auto dimensionality(T const(&t)[N]) {return 1 + dimensionality(t[0]);}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T, typename = decltype(std::declval<T>().sizes())>
       auto has_sizes_aux(T const&) -> std::true_type;
inline auto has_sizes_aux(...     ) -> std::false_type;

template<class T> struct has_sizes : decltype(has_sizes_aux(std::declval<T>())){};

template<class Array, typename = std::enable_if_t<has_sizes<Array>{}> >
constexpr auto sizes(Array const& arr)
->decltype(arr.sizes()) {
	return arr.sizes(); }

template<class T, typename = std::enable_if_t<not has_sizes<T>{}> >
inline constexpr auto sizes(T const& /*unused*/) -> std::tuple<> {return {};}

template<class T, std::size_t N>
constexpr auto sizes(const T(&t)[N]) noexcept {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
//  using std::size; // this line needs c++17
	using boost::multi::size;
	return std::tuple_cat(std::make_tuple(boost::multi::size(t)), sizes(t[0]));
}

template<class T, std::size_t N>
constexpr auto base(T(&t)[N]) noexcept {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	return data_elements(t);
}

template<class T, std::size_t N>
constexpr auto base(T(*&t)[N]) noexcept {return base(*t);}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T, typename = std::enable_if_t<not std::is_array<T>{}> >
constexpr auto base(T const* t) noexcept {return t;}

template<class T, typename = std::enable_if_t<not std::is_array<T>{}> >
constexpr auto base(T* t) noexcept {return t;}

template<class T, std::enable_if_t<std::is_standard_layout<T>{} and std::is_trivial<T>{}, int> = 0>
auto base(T& t) {return &t;}

template<class T>
constexpr auto corigin(const T& t) {return &t;}

template<class T, std::size_t N>
constexpr auto corigin(const T(&t)[N]) noexcept {return corigin(t[0]);}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility

template<class T, typename = decltype(std::declval<T>().extension())>
       auto has_extension_aux(T const&) -> std::true_type;
inline auto has_extension_aux(...     ) -> std::false_type;
template<class T> struct has_extension : decltype(has_extension_aux(std::declval<T>())){};

template<class Container, class=std::enable_if_t<not has_extension<Container>{}>>
auto extension(Container const& c)  // TODO(correaa) consider "extent"
->decltype(multi::extension_t<std::make_signed_t<decltype(size(c))>>(0, static_cast<std::make_signed_t<decltype(size(c))>>(size(c)))) {
	return multi::extension_t<std::make_signed_t<decltype(size(c))>>(0, static_cast<std::make_signed_t<decltype(size(c))>>(size(c))); }

template<class T, typename = decltype(std::declval<T>().extensions())>
       auto has_extensions_aux(T const&) -> std::true_type;
inline auto has_extensions_aux(...     ) -> std::false_type;

template<class T> struct has_extensions : decltype(has_extensions_aux(std::declval<T>())) {};

template<class T, typename = std::enable_if_t<has_extensions<T>{}> >
NODISCARD("") auto extensions(T const& t)
->std::decay_t<decltype(t.extensions())> {
	return t.extensions();
}

template<class T, typename = std::enable_if_t<not has_extensions<T>{}> >
constexpr auto extensions(T const& /*unused*/) -> multi::layout_t<0>::extensions_type {return {};}

template<class T, std::size_t N>
constexpr auto extensions(T(&t)[N]) {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : for backwards compatibility
	return index_extension(N)*extensions(t[0]);
}

template<dimensionality_type D>
struct extensions_aux {
	template<class T>
	static auto call(T const& t) {
		return std::tuple_cat(std::make_tuple(t.extension()), extensions<D-1>(t));
	}
};

template<> struct extensions_aux<0> {
	template<class T> static auto call(T const& /*unused*/){return std::make_tuple();}
};

template<dimensionality_type D, class T>
auto extensions(T const& t) {
	return extensions_aux<D>::call(t);
}

template<class T1> struct extensions_t_aux;

template<class T1, class T2> auto extensions_(T2 const& t2) {
	return extensions_t_aux<T1>::call(t2);
}

template<class T1> struct extension_t_aux {
	static auto call(T1 const& /*unused*/) {return std::make_tuple();}
	template<class T2>
	static auto call(T2 const& t2) {return std::tuple_cat(std::make_tuple(t2.extension()), extensions_<T1>(*begin(t2)));}
};

template<class T, typename = decltype(std::declval<T const&>().layout())>
       auto has_layout_member_aux(T const&) -> std::true_type;
inline auto has_layout_member_aux(...     ) -> std::false_type;

template<class T>
struct has_layout_member : decltype(has_layout_member_aux(std::declval<T const&>())){};

template<class T, typename = std::enable_if_t<has_layout_member<T const&>{}> >
auto layout(T const& t)
->decltype(t.layout()) {
	return t.layout(); }

template<class T, typename = std::enable_if_t<not has_layout_member<T const&>{}> >
auto layout(T const& /*unused*/) -> layout_t<0> {return {};}

template<class T, std::size_t N>
constexpr auto layout(T(&t)[N]) {return multi::layout_t<std::rank<T[N]>{}>(multi::extensions(t));}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): for backward compatibility

template<class T, std::size_t N>
constexpr auto strides(T(&t)[N]) {return layout(t).strides();}  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): for backward compatibility

template<class T, std::size_t N>
struct array_traits<std::array<T, N>> {
	static constexpr auto dimensionality() -> dimensionality_type {return 1;}

	using reference = T&;
	using value_type = std::decay_t<T>;
	using pointer = T*;
	using element = value_type;
	using element_ptr = pointer;
	using decay_type = std::array<value_type, N>;
};

template<class T, std::size_t N, std::size_t M>
struct array_traits<std::array<std::array<T, M>, N>> {
	static constexpr auto dimensionality() -> dimensionality_type {return 1 + array_traits<std::array<T, M>>::dimensionality();}

	using reference = std::array<T, M>&;
	using value_type = std::array<std::decay_t<T>, M>;
	using pointer = std::array<T, M>*;
	using element = typename array_traits<std::array<T, M>>::element;
	using element_ptr = typename array_traits<std::array<T, M>>::element;
	using decay_type = std::array<value_type, M>;
};

template<class T, std::size_t N> constexpr auto data_elements(std::array<T, N>&       arr) noexcept {return arr.data();}
template<class T, std::size_t M, std::size_t N> constexpr auto data_elements(std::array<std::array<T, M>, N>& arr) noexcept {return data_elements(arr[0]);}

template<class T, std::size_t N> constexpr auto data_elements(std::array<T, N> const&       arr) noexcept {return arr.data();}
template<class T, std::size_t M, std::size_t N> constexpr auto data_elements(std::array<std::array<T, M>, N> const& arr) noexcept {return data_elements(arr[0]);}

template<class T, std::size_t N> constexpr auto data_elements(std::array<T, N>&&       arr) noexcept {return arr.data();}

template<class T, std::size_t M, std::size_t N>
constexpr auto data_elements(std::array<std::array<T, M>, N>&& arr) noexcept {return data_elements(arr[0]);}

template<class T, std::size_t N> constexpr auto num_elements(std::array<T, N> const& /*unused*/) noexcept
-> std::ptrdiff_t{return N;}

template<class T, std::size_t M, std::size_t N>
constexpr auto num_elements(std::array<std::array<T, M>, N> const& arr)
-> std::ptrdiff_t {return static_cast<std::ptrdiff_t>(N)*num_elements(arr[0]);}

template<class T, std::size_t N>
constexpr auto dimensionality(std::array<T, N> const& /*unused*/) -> dimensionality_type {return 1;}

template<class T, std::size_t M, std::size_t N>
constexpr auto dimensionality(std::array<std::array<T, M>, N> const& arr) -> dimensionality_type {
	return 1 + dimensionality(arr[0]);
}

#if (__cplusplus < 201703L)
// this conflicts with std::size in nvcc 11 and c++17
template<class T, std::size_t N>
constexpr auto size(std::array<T, N> const& /*arr*/) {
	return multi::size_type{N};
}
#endif

template<class T, std::size_t N>
constexpr auto extensions(std::array<T, N> const& /*arr*/) {
	return multi::extensions_t<1>{{0, N}};
}

template<class T, std::size_t N, std::size_t M>
auto extensions(std::array<std::array<T, N>, M> const& arr) {
	return multi::iextension(M)*extensions(arr[0]);
}

template<class T, std::size_t N>
constexpr auto stride(std::array<T, N> const& /*arr*/) {
	return multi::size_type{1};  // multi::stride_type?
}

template<class T, std::size_t N, std::size_t M>
constexpr auto stride(std::array<std::array<T, N>, M> const& arr) {
	return num_elements(arr[0]);
}

template<class T, std::size_t N>
constexpr auto layout(std::array<T, N> const& arr) {
	return multi::layout_t<multi::array_traits<std::array<T, N>>::dimensionality()>{multi::extensions(arr)};
}

}  // end namespace multi
}  // end namespace boost

#endif

