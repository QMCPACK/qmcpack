// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa

//#ifndef MULTI_CONFIG_MARK_HPP
//#define MULTI_CONFIG_MARK_HPP

#ifndef MULTI_MARK_SCOPE  // NOLINT(llvm-header-guard) this is a configuration header, can be included many times
	#ifdef CALI_CXX_MARK_SCOPE
		#define MULTI_MARK_SCOPE(MsG) CALI_CXX_MARK_SCOPE(MsG)
	#else
		#define MULTI_MARK_SCOPE(MsG) ((void)0)  // NOLINT(cppcoreguidelines-macro-usage) to mark scopes
	#endif
#endif

//#endif
