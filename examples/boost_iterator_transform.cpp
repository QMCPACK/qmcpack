#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $CXXFLAGS $0 -o $0.$X&&$0.$X&&rm $0.$X;exit
#endif
//  © Alfredo A. Correa 2018-2020

#include "../array.hpp"

#include<boost/iterator/transform_iterator.hpp>
#include<boost/functional/hash.hpp>

namespace multi = boost::multi;

int main(){

{
	auto r = multi::make_range(5, 10);
	auto f = [](auto x){return x+1;};
	std::vector<double> v(
		boost::make_transform_iterator(r.begin(), f), 
		boost::make_transform_iterator(r.end()  , f)
	);
	assert( v[1] == 7 );
}
{
	auto r = multi::make_range(5, 10);
	auto f = [](auto x){return x+1;};
	multi::array<double, 1> v(
		boost::make_transform_iterator(r.begin(), f), 
		boost::make_transform_iterator(r.end()  , f)
	);
	assert( v[1] == 7 );
}
{
	multi::array<double, 1> v(10);
	auto r = extension(v);
	auto f = [](auto x){return x*2;};
	v.assign(
		boost::make_transform_iterator(r.begin(), f),
		boost::make_transform_iterator(r.end()  , f)
	);
	assert( v[1] == 2 );
}
{
	auto r = multi::make_extension_t(10l);
	auto f = [](auto x){
		std::size_t seed = 1234;
	//	boost::hash_combine(seed, );
		seed ^= boost::hash<multi::index>{}(x) + 0x9e3779b9 + (seed<<6) + (seed>>2);
		return static_cast<double>(seed)/static_cast<double>(std::numeric_limits<std::size_t>::max());
	};
	multi::array<double, 1> v(
		boost::make_transform_iterator(r.begin(), f), 
		boost::make_transform_iterator(r.end()  , f)
	);

	std::size_t seed = 12349l;
	//	boost::hash_combine(seed, );
//	seed ^= boost::hash<std::size_t>{}(13) + 0x9e3779b9 + (seed<<6) + (seed>>2);
	boost::hash_combine(seed, 13);

	assert( v.size() == r.size() );
	assert( v[1] >= 0. );
	assert( v[1] < 1.  );
	assert( std::all_of(begin(v), end(v), [](auto x){
		return x >= 0. and x < 1.;
	}) );
}

}

