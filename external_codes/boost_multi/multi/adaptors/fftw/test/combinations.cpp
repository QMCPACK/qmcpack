#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lfftw3 -lboost_unit_test_framework -lboost_timer&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW adaptor (cpu) with thrust complex"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include <boost/timer/timer.hpp>

#include "../../fftw.hpp"

#include<complex>
#include<chrono>

#include<thrust/complex.h>

namespace multi = boost::multi;

namespace utf = boost::unit_test::framework;

struct watch : private std::chrono::high_resolution_clock{
	std::string label_; time_point  start_;
	watch(std::string label ="") : label_{label}, start_{now()}{}
	~watch(){
		std::cerr<< label_<<": "<< std::chrono::duration<double>(now() - start_).count() <<" sec"<<std::endl;
	}
};

BOOST_AUTO_TEST_CASE(fft_combinations, *boost::unit_test::tolerance(0.00001)){

	using complex = std::complex<double>;

	auto const in = []{
		multi::array<complex, 4> ret({32, 90, 98, 96});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();
	std::cout<<"memory size "<< in.num_elements()*sizeof(complex)/1e6 <<" MB\n";

	std::vector<std::array<bool, 4>> cases = {
		{false, true , true , true }, 
		{false, true , true , false}, 
		{true , false, false, false}, 
		{true , true , false, false},
		{false, false, true , false},
		{false, false, false, false},
	};

	using std::cout;
	for(auto c : cases){
		cout<<"case "; copy(begin(c), end(c), std::ostream_iterator<bool>{cout,", "}); cout<<"\n";
		multi::array<complex, 4> out = in;
		{
			boost::timer::auto_cpu_timer t{"cpu_oplac %ws wall, CPU (%p%)\n"};
			multi::fftw::dft_forward(c, in, out);
		}
		{
			multi::fftw::plan p(c, in, out, multi::fftw::forward);
			boost::timer::auto_cpu_timer t{"cpu_oplac planned %ws wall, CPU (%p%)\n"};
			p();
		}
		{
			auto in_rw = in;
			boost::timer::auto_cpu_timer t{"cpu_iplac %ws wall, CPU (%p%)\n"};
			multi::fftw::dft_forward(c, in_rw);
		//	BOOST_TEST( abs( in_rw[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
		{
			auto in_rw = in;
			multi::fftw::plan p(c, in_rw, in_rw, multi::fftw::forward);
			boost::timer::auto_cpu_timer t{"cpu_iplac planned %ws wall, CPU (%p%)\n"};
			p();
		//	BOOST_TEST( abs( in_rw[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
		{
			auto in_rw = in;
			multi::fftw::plan p(c, in_rw, in_rw, multi::fftw::forward);// | FFTW_MEASURE);
			boost::timer::auto_cpu_timer t{"cpu_iplac planned measured %ws wall, CPU (%p%)\n"};
			p();
		//	BOOST_TEST( abs( in_rw[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
		{
			boost::timer::auto_cpu_timer t{"cpu_alloc %ws wall, CPU (%p%)\n"}; 
			auto out_cpy = multi::fftw::dft_forward(c, in);
			BOOST_TEST( abs( out_cpy[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
		{
			auto in_rw = in;
			boost::timer::auto_cpu_timer t{"cpu_move %ws wall, CPU (%p%)\n"}; 
			auto out_cpy = multi::fftw::dft_forward(c, std::move(in_rw));
			BOOST_REQUIRE( in_rw.empty() );
			BOOST_TEST( abs( out_cpy[5][4][3][1] - out[5][4][3][1] ) == 0. );
		}
	}
}

BOOST_AUTO_TEST_CASE(fftw_4D_power_benchmark, *boost::unit_test::disabled() ){
	using complex = std::complex<double>;
	namespace fftw = multi::fftw;

	auto x = multi::array<complex, 4>::extensions_type({64, 128, 128, 128});
	multi::array<complex, 4> in(x);
	std::iota(in.data_elements(), in.data_elements() + in.num_elements(), 1.2);

	BOOST_REQUIRE( in[0][0][0][0] == 1.2 );
	std::array<bool, 4> c = {false, true, true, true};
	[&, _ = watch{utf::current_test_case().full_name()+" inplace FTTT"}]{
		fftw::dft(c, in, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" inplace FTTT"}]{
		fftw::dft(c, in, fftw::forward);
	}();
	auto in0000 = in[0][0][0][0];
	BOOST_REQUIRE( in0000 != 1.2 );


	multi::array<complex, 4> out(x);
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace FTTT"}]{
		fftw::dft(c, in, out, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace FTTT"}]{
		fftw::dft(c, in, out, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace FTTT"}]{
		fftw::dft(c, in, out, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace+alloc FTTT"}]{
		multi::array<complex, 4> out2(x);
		fftw::dft(c, in, out2, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace+alloc FTTT"}]{
		multi::array<complex, 4> out2(x);
		fftw::dft(c, in, out2, fftw::forward);
	}();
	BOOST_REQUIRE( in0000 == in[0][0][0][0] );
	
}

