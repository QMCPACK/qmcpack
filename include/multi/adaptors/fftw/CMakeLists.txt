cmake_minimum_required(VERSION 3.11)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#project(
#	boost-multi-adaptors-fftw
#	VERSION 0.1
#	LANGUAGES CXX
#)

find_package(PkgConfig)

if(PKG_CONFIG_FOUND)
	pkg_search_module(
		FFTW
		fftw3
		IMPORTED_TARGET
	)
	if(FFTW_FOUND)
		include_directories(PkgConfig::FFTW)
		link_libraries(PkgConfig::FFTW)

		include_directories(${CMAKE_BINARY_DIR})

		add_subdirectory(test)
	else()
		message(WARNING "Cannot find FFTW, FFTW-adaptor will not be tested. If you want this feature install FFTW, for example please run:\n    sudo apt install pkg-config libfftw3-dev")
	endif()
else()
	message(WARNING "Cannot find PkgConfig and/or FFTW, FFTW-adaptor will not be tested. If you want this feature install PkgConfig and FFTW, for example please run:\n    sudo apt install pkg-config libfftw3-dev")
endif()

