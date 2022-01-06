#ifdef COMPILATION_INSTRUCTIONS//-*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*-
$CXXX $CXXFLAGS $0 -o $0x `pkg-config --cflags --libs cudart-11.1`&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_MEMORY_ADAPTOR_CUDA_DETAIL_ERROR_HPP
#define MULTI_MEMORY_ADAPTOR_CUDA_DETAIL_ERROR_HPP

#include<driver_types.h>     // cudaError_t
#include<cuda_runtime_api.h> // cudaGetErrorString

#include<system_error>
#include<type_traits>        // underlying_type

namespace Cuda{

enum /*class*/ error : std::underlying_type<cudaError_t>::type{
	success                    = cudaSuccess, // = 0 The API call returned with no errors. In the case of query calls, this also means that the operation being queried is complete (see cudaEventQuery() and cudaStreamQuery()).
	missing_configuration      = cudaErrorMissingConfiguration,
//	invalid_value /*invalid_argument*/ = cudaErrorInvalidValue, // = 1, This indicates that one or more of the parameters passed to the API call is not within an acceptable range of values.
	memory_allocation          = cudaErrorMemoryAllocation, // = 2 // The API call failed because it was unable to allocate enough memory to perform the requested operation. 
	initialization_error       = cudaErrorInitializationError, 
	lauch_failure              = cudaErrorLaunchFailure,
	lauch_timeout              = cudaErrorLaunchTimeout, 
	lauch_out_of_resources     = cudaErrorLaunchOutOfResources, 
	invalid_device_function    = cudaErrorInvalidDeviceFunction, 
	invalid_configuration      = cudaErrorInvalidConfiguration, 
	invalid_device             = cudaErrorInvalidDevice, 
	invalid_value              = cudaErrorInvalidValue, ///*invalid_argument*/ = cudaErrorInvalidValue, // = 1 This indicates that one or more of the parameters passed to the API call is not within an acceptable range of values. 
	invalid_pitch_value        = cudaErrorInvalidPitchValue, 
	invalid_symbol             = cudaErrorInvalidSymbol, 
	unmap_buffer_object_failed = cudaErrorUnmapBufferObjectFailed, 
	invalid_device_pointer     = cudaErrorInvalidDevicePointer, 
	invalid_texture            = cudaErrorInvalidTexture, 
	invalid_texture_binding    = cudaErrorInvalidTextureBinding, 
	invalid_channel_descriptor = cudaErrorInvalidChannelDescriptor, 
	invalid_memcpy_direction   = cudaErrorInvalidMemcpyDirection, 
	invalud_filter_setting     = cudaErrorInvalidFilterSetting, 
	invalid_norm_setting       = cudaErrorInvalidNormSetting, 
	unknown                    = cudaErrorUnknown, 
	invalid_resource_handle    = cudaErrorInvalidResourceHandle, 
	insuffient_driver          = cudaErrorInsufficientDriver, 
	no_device                  = cudaErrorNoDevice, 
	set_on_active_process      = cudaErrorSetOnActiveProcess, 
	startup_failure            = cudaErrorStartupFailure, 
	invalid_ptx                = cudaErrorInvalidPtx, 
	no_kernel_image_for_device = cudaErrorNoKernelImageForDevice, 
	jit_compiler_not_found     = cudaErrorJitCompilerNotFound
};

inline std::string string(enum error e){return cudaGetErrorString(static_cast<cudaError_t>(e));}

struct error_category : std::error_category{
	char const* name() const noexcept override{return "cuda wrapper";}
	std::string message(int e) const override{return string(static_cast<error>(e));}
	static error_category& instance(){
		static error_category instance;
		return instance;
	}
};

inline std::error_code make_error_code(error err) noexcept{
	return {int(err), error_category::instance()};
}

}

namespace std{template<> struct is_error_code_enum<Cuda::error> : true_type{};}

#if not __INCLUDE_LEVEL__

#include<iostream>

using std::cout;

int main(){

	{
		std::error_code ec = Cuda::error::memory_allocation; (void)ec;
	}
	try{
		auto e = Cuda::error::memory_allocation; // return from a cudaFunction
		throw std::system_error{e, "I cannot do allocation"};
	}catch(std::system_error const& e){
		cout
			<<"catched...\n"
			<<"code: "   << e.code()           <<'\n'
			<<"message: "<< e.code().message() <<'\n'
			<<"what: "   << e.what()           <<'\n'
		;
	}

//	auto e = Cuda::error::memory_allocation; // return from a cudaFunction
//	throw std::system_error{e, "because"};

}
#endif
#endif

