#pragma once

#include<driver_types.h>     // cudaError_t
#include<cuda_runtime_api.h> // cudaGetErrorString

#include<system_error>
#include<type_traits>        // underlying_type

namespace boost{
namespace multi{
namespace cuda{
namespace runtime{

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
	return {static_cast<int>(err), error_category::instance()};
}

}}}}

namespace std{template<> struct is_error_code_enum<boost::multi::cuda::runtime::error> : true_type{};}

