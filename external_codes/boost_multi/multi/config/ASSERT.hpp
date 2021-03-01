#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXXX $CXXFLAGS -DMULTI_NDEBUG_ACCESS $0 -o $0x &&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_CONFIG_ASSERT_HPP
#define MULTI_CONFIG_ASSERT_HPP

#include<cassert>

#ifdef MULTI_ACCESS_NDEBUG
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


