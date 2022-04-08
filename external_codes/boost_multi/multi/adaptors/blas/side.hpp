#ifndef MULTI_ADAPTORS_BLAS_SIDE_HPP // -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
#define MULTI_ADAPTORS_BLAS_SIDE_HPP
// © Alfredo A. Correa 2019-2021

namespace boost{
namespace multi{
namespace blas{

//enum class SIDE : char{L='L', R='R'};

enum side : char{
	left  = 'L', 
	right = 'R'//,
//	pre_multiply = 'R', 
//	post_multiply = 'L'
};

inline auto swap(side s) -> side{
	switch(s){
		case side::left: return side::right;
		case side::right: return side::left;
	} __builtin_unreachable();
}

} // end namespace blas
} // end namespace multi
} // end namespace boost

#endif

