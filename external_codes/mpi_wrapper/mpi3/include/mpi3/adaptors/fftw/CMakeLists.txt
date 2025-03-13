if(PKG_CONFIG_FOUND)
	pkg_search_module(
		FFTW
		fftw3
		IMPORTED_TARGET
	)
	if(FFTW_FOUND)
        add_subdirectory(test)
	else()
		message(WARNING "Cannot find FFTW, FFTW-adaptor will not be tested. If you want this feature install MPI FFTW, for example please run:"
		                "\n    sudo apt install pkg-config libfftw3-dev libfftw3-mpi-dev"
		                "\n    sudo dnf install fftw-devel  # in Fedora")
	endif()
else()
    message(WARNING "Cannot find FFTW, FFTW-adaptor will not be tested. If you want this feature install MPI FFTW, for example please run:"
            "\n    sudo apt install pkg-config libfftw3-dev libfftw3-mpi-dev"
            "\n    sudo dnf install fftw-devel  # in Fedora")
endif()
