// Copyright 2020-2024 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_HIPFFT_HPP
#define MULTI_ADAPTORS_HIPFFT_HPP

#include <hipfft/hipfft.h>
#include <hipfft/hipfftXt.h>


using cudaError_t = hipError_t;

constexpr static auto const& cudaDeviceReset  = hipDeviceReset;
constexpr static auto const& cudaDeviceSynchronize  = hipDeviceSynchronize;
constexpr static auto const& cudaSuccess = hipSuccess;

#define cu2hip_fft(TypeleafnamE) using cufft ## TypeleafnamE = hipfft ## TypeleafnamE
    cu2hip_fft(Handle);
    cu2hip_fft(DoubleComplex);
    cu2hip_fft(Result);
#undef cu2hip_fft

#define cu2hip_fft(FunctionleafnamE) constexpr static auto const& cufft ## FunctionleafnamE  = hipfft ## FunctionleafnamE
    cu2hip_fft(Create);
    cu2hip_fft(Destroy);
    cu2hip_fft(GetSize);
    cu2hip_fft(ExecZ2Z);
    cu2hip_fft(SetAutoAllocation);
    cu2hip_fft(SetWorkArea);
    cu2hip_fft(PlanMany);
#undef cu2hip_fft

#define CU2HIPFFT_(NamE) constexpr static auto const& CUFFT_ ## NamE  = HIPFFT_ ## NamE

CU2HIPFFT_(ALLOC_FAILED);
CU2HIPFFT_(BACKWARD);

constexpr static auto const& CUFFT_INVERSE = HIPFFT_BACKWARD;

CU2HIPFFT_(EXEC_FAILED);
CU2HIPFFT_(FORWARD);
CU2HIPFFT_(INCOMPLETE_PARAMETER_LIST);
CU2HIPFFT_(INTERNAL_ERROR);
CU2HIPFFT_(INVALID_DEVICE);
CU2HIPFFT_(INVALID_SIZE);
CU2HIPFFT_(INVALID_TYPE);
CU2HIPFFT_(INVALID_VALUE);
CU2HIPFFT_(INVALID_PLAN);
CU2HIPFFT_(NO_WORKSPACE);
CU2HIPFFT_(NOT_IMPLEMENTED);
CU2HIPFFT_(NOT_SUPPORTED);
CU2HIPFFT_(UNALIGNED_DATA);
CU2HIPFFT_(PARSE_ERROR);
CU2HIPFFT_(SETUP_FAILED);
CU2HIPFFT_(SUCCESS);
CU2HIPFFT_(Z2Z);

#undef CU2HIPFFT_

#include "cufft.hpp"

// namespace boost::multi{
//     namespace cufft = hipfft;
// }

#endif
