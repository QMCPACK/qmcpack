// Copyright 2022-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

// #include <boost/multi/adaptors/fftw.hpp>  // includes fftw3.hpp

#include <complex>
#include <functional>
#include <numeric>  // for std::iota

namespace multi = boost::multi;

auto main() -> int {
	using complex = std::complex<double>;

	// input array
	auto const x = std::invoke([] {  // NOLINT(readability-identifier-length)
		multi::array<complex, 1> ret(8);
		// fill the first array with some numbers
		std::iota(ret.begin(), ret.end(), 1.0);
		return ret;
	});

	// output array
	// multi::array<complex, 1> y(x.size());  // NOLINT(readability-identifier-length)
	// compute the FFT of x and store results in y
	// auto y = +multi::fftw::dft_forward(x);  // NOLINT(readability-identifier-length)

	// display the results
	// std::cout << "FFT =" << std::endl;
	// std::copy(y.begin(), y.end(), std::ostream_iterator<complex>(std::cout, "\n"));

	// "shifted" results
	// std::rotate(y.begin(), y.begin() + y.size() / 2 + y.size() % 2, y.end());

	// std::cout << "FFT shifted =" << std::endl;
	// std::copy(y.begin(), y.end(), std::ostream_iterator<complex>(std::cout, "\n"));
}
