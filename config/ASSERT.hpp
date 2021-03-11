// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
// © Alfredo A. Correa 2019-2021

#ifndef MULTI_CONFIG_ASSERT_HPP
#define MULTI_CONFIG_ASSERT_HPP

#include<cassert>

#if defined(MULTI_ACCESS_NDEBUG) or defined(__CUDACC__)
#define MULTI_ACCESS_ASSERT(Expr)
#else
#define MULTI_ACCESS_ASSERT(Expr) assert(Expr)
#endif

#if not __INCLUDE_LEVEL__

int main(){
	MULTI_ACCESS_ASSERT(false && "hola");
}
#endif
#endif

