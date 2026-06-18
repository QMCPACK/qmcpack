#if COMPILE_RUN_INSTRUCTIONS
${CXX:-c++} -std=c++2b $0 -I../include && ./a.out; exit
#endif
#include <boost/multi/array.hpp>

#include <cmath>  // for std::exp
#include <iostream>
#include <ranges>

namespace multi = boost::multi;
namespace stv   = std::views;

void plot(auto const& x, auto const& f, std::string const& title = "") {
	assert(x.size() == f.size());
	std::cout << "set title '" << title << "'\n"
			  << "plot '-' with linespoints\n";
	for(auto i : x.extension()) {
		std::cout << x[i] << " " << f[i] << "\n";
	}
	std::cout << 'e' << std::endl
			  << "pause 0.1\n";
}

template<class Range, typename T>
auto append(Range const& range, T const& value) {
	return stv::iota(typename Range::size_type{}, range.size() + 1) | stv::transform([&](auto i) -> decltype(auto) {return (i<range.size()-1)?range[i]:value;});
}

template<class Range, typename T>
auto prepend(Range const& range, T const& value) {
	return stv::iota(typename Range::size_type{}, range.size() + 1) | stv::transform([&](auto i) -> decltype(auto) {return (i==0)?value:range[i-1];});
}

auto main() -> int {

	using multi::operator+;

	// dx = 0.2;
	// x = [0:1:20]*dx;
	// f = x.*exp(-x.^2);
	// plot(x, f, 'ro');

	auto dx = 0.2;
	auto x  = +(stv::iota(0, 20) | stv::transform([dx](auto i) { return i * dx; }));
	auto f  = +(x | stv::transform([](auto e) { return e * std::exp(-e * e); }));
	plot(x, f);

	// f_my_left = [NaN, f(1:end-1)];
	// f_my_right = [f(2:end), NaN];
	// d2f = (f_my_right - 2*f + f_my_left)/(dx^2);
	auto f_my_left = +prepend(f.taked(f.size() - 1), NAN);
	auto f_my_right = +append(f.dropped(1), NAN);
	auto d2f = +stv::zip_transform([dx2 = dx * dx](auto r, auto m, auto l) { return (r - 2 * m + l) / dx2; }, f_my_right, f, f_my_left);

	// dt = 0.01; D = 1;
	// for k=1:100,
	//     f_my_left = [NaN, f(1:end-1)];
	//     f_my_right = [f(2:end), NaN];
	//     d2f = (f_my_right - 2*f + f_my_left)/(dx^2);
	//     f(2:end-1) = f(2:end-1) + D*dt*d2f(2:end-1);
	//     plot(x, f, 'ro-'); ylim([0 0.45]); drawnow
	//     pause(0.1)
	// end

	auto dt = 0.01;
	auto D  = 1.0;
	for(auto k = 0; k != 100; ++k) {
		f_my_left({1, f.size()})      = f({0, f.size() - 1});
		f_my_right({0, f.size() - 1}) = f({1, f.size()});

		d2f                  = stv::zip_transform([dx2 = dx * dx](auto r, auto m, auto l) { return (r - 2 * m + l) / dx2; }, f_my_right, f, f_my_left);
		f({1, f.size() - 1}) = stv::zip_transform([&](auto eff, auto d2) { return eff + D * dt * d2; }, f({1, f.size() - 1}), d2f({1, f.size() - 1}));
		plot(x, f, "k=" + std::to_string(k));
	}
}
