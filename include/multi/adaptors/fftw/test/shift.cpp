// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW shift"
#include<boost/test/unit_test.hpp>

#include "../../fftw.hpp"
#include <random>

template<class T>
class n_random_complex {
	std::size_t n_ = 0;
    mutable std::mt19937 gen_{std::random_device{}()};
	mutable std::uniform_real_distribution<> dist_{-1., 1.};

 public:
	explicit n_random_complex(std::size_t n) : n_{n} {}

	class iterator : public boost::multi::random_access_iterator<iterator, std::complex<T>, std::complex<T>, void> {
		n_random_complex<T> const* ptr_;
		std::size_t n_;

	 public:
		iterator(n_random_complex<T> const* ptr, std::size_t n) : ptr_{ptr}, n_{n} {}

		auto operator*() const {return std::complex<T>{ptr_->dist_(ptr_->gen_), ptr_->dist_(ptr_->gen_)};}
		auto operator++() -> iterator& {++n_; return *this;}
		friend auto operator==(iterator const& s, iterator const& o) {return s.n_ == o.n_;}
		auto operator-(iterator const& other) const {return n_ - other.n_;}
	};
	auto begin() const {return iterator{this, 0 };}
	auto end  () const {return iterator{this, n_};}
	auto size()  const {return n_;}
};

namespace multi = boost::multi;
namespace fftw = multi::fftw;

BOOST_AUTO_TEST_CASE(fftw_shift){

	class watch : std::chrono::steady_clock {
		time_point start_ = now();
		public: auto elapsed_sec() const {return std::chrono::duration<double>(now() - start_).count();}
	};

	multi::array<std::complex<double>, 1> const arr = n_random_complex<double>(19586);  BOOST_REQUIRE(arr.size() == 19586);
	multi::array<std::complex<double>, 1>       res(arr.extensions());                  BOOST_REQUIRE(res.size() == 19586);

	fftw::plan fdft{arr, res, multi::fftw::forward};

	auto const N = 40;
	[&, _ = watch{}]{
		for(int i = 0; i != N; ++i) {
			fdft(arr.base(), res.base());
			std::rotate(res.begin(), res.begin() + res.size()/2, res.end());
		}

	    BOOST_TEST_MESSAGE( "FFTW shift " << _.elapsed_sec()/N <<" sec" );  // prints  0.000882224 sec
	}();
}
