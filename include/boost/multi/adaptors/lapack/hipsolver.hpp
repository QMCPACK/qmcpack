// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_LAPACK_HIPSOLVER_HPP
#define BOOST_MULTI_ADAPTORS_LAPACK_HIPSOLVER_HPP
#pragma once

// include this header first (or compile with -DMULTI_USE_HIP globally, as the ROCm CI does)
// so that the other Multi adaptors (thrust, cublas) also take their HIP branches

#if !defined(MULTI_USE_HIP)
#define MULTI_USE_HIP
#endif

#include "cusolver.hpp"  // under MULTI_USE_HIP this is the hipSOLVER backend

#endif  // BOOST_MULTI_ADAPTORS_LAPACK_HIPSOLVER_HPP
