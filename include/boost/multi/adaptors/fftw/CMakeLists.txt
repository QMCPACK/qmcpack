cmake_minimum_required(VERSION 3.11)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

add_library(multi_fftw3 INTERFACE)

find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
	pkg_search_module(
		FFTW
		fftw3
		IMPORTED_TARGET
	)
	if(FFTW_FOUND)
		# include_directories(PkgConfig::FFTW) this seems to be always incorrect
		link_libraries(PkgConfig::FFTW)
		target_link_libraries(multi_fftw3 INTERFACE PkgConfig::FFTW)

		include_directories(${CMAKE_BINARY_DIR})

		add_subdirectory(test)
#       add_subdirectory(mpi)
	else()
		message(WARNING "Cannot find FFTW, FFTW-adaptor will not be tested. If you want this feature install FFTW, for example please run:"
		                "\n    sudo apt install pkg-config libfftw3-dev"
		                "\n    sudo dnf install fftw-devel  # in Fedora"
						"\n    brew install fftw pkgconf")
	endif()
else()
	message(WARNING "Cannot find PkgConfig and/or FFTW, FFTW-adaptor will not be tested. If you want this feature install PkgConfig and FFTW, for example please run:"
					"\n    sudo apt install pkg-config libfftw3-dev"
					"\n    sudo dnf install fftw-devel  # in Fedora")

	find_package(FFTW3 COMPONENTS fftw3)  # or find_package(FFTW3 COMPONENTS fftw3f fftw3) for specific components

	if(FFTW3_FOUND)
		message(WARNING "Found FFTW3")
		target_link_libraries(multi_fftw3 INTERFACE FFTW3::fftw3)
		add_subdirectory(test)

	else()
		message(WARNING "Cannot find FFTW3, FFTW3-adaptor will not be tested. If you want this feature install PkgConfig and FFTW, for example please run:"
						"\n    sudo apt install libfftw3-dev"
						"\n    sudo dnf install fftw-devel  # in Fedora"
						"\n    vcpkg install fftw3  # on Windows")
	endif()
endif()
