find_package(PkgConfig)

set(CMAKE_CXX_EXTENSIONS OFF)

if(PKG_CONFIG_FOUND)
	pkg_search_module(
		FFTW
		fftw3
		IMPORTED_TARGET
	)
	if(FFTW_FOUND)
		# include_directories(PkgConfig::FFTW) this seems to be always incorrect
		link_libraries(PkgConfig::FFTW)

		include_directories(${CMAKE_BINARY_DIR})

		#add_subdirectory(test)
	else()
		message(WARNING "Cannot find FFTW, FFTW-adaptor will not be tested. If you want this feature install FFTW, for example please run:"
		                "\n    sudo apt install pkg-config libfftw3-dev"
		                "\n    sudo dnf install fftw-devel  # in Fedora")
	endif()
else()
	message(WARNING "Cannot find PkgConfig and/or FFTW, FFTW-adaptor will not be tested. If you want this feature install PkgConfig and FFTW, for example please run:"
					"\n    sudo apt install pkg-config libfftw3-dev"
					"\n    sudo dnf install fftw-devel  # in Fedora")
endif()

include(FetchContent)
FetchContent_Declare(multi GIT_REPOSITORY https://gitlab.com/correaa/boost-multi.git)
FetchContent_MakeAvailable(multi)

set(TEST_SRCS
  array_2d.cpp
)
  
foreach(TEST_FILE ${TEST_SRCS})
  set(TEST_EXE "${TEST_FILE}.x")
  add_executable(${TEST_EXE} ${TEST_FILE})
  target_link_libraries(${TEST_EXE} PRIVATE bmpi3 multi PkgConfig::FFTW fftw3_mpi)
  add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
endforeach()
