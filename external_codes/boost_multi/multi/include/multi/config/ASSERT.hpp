// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#ifndef MULTI_CONFIG_ASSERT_HPP
#define MULTI_CONFIG_ASSERT_HPP

#include<cassert>

#if defined(MULTI_ACCESS_NDEBUG) or defined(__CUDACC__)
	#define MULTI_ACCESS_ASSERT(Expr)
#else
	// NOLINTNEXTLINE(cppcoreguidelines-macro-usage) this is for very inefficient asserts
	#define MULTI_ACCESS_ASSERT(Expr) assert(Expr)
#endif

#endif
