// Copyright 2020-2024 Alfredo A. Correa

#include <fstream>

#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/nvp.hpp>

#include <boost/archive/polymorphic_binary_iarchive.hpp>
#include <boost/archive/polymorphic_binary_oarchive.hpp>
#include <boost/archive/polymorphic_text_iarchive.hpp>
#include <boost/archive/polymorphic_text_oarchive.hpp>
#include <boost/archive/polymorphic_xml_iarchive.hpp>
#include <boost/archive/polymorphic_xml_oarchive.hpp>

#include "../../multi/array.hpp"
#include <experimental/filesystem>

enum format { xml,
              txt,
              bin };

namespace barch = boost::archive;
namespace bs11n = boost::serialization;

#define BOOST_MULTI_UNSWITCH __builtin_unreachable

template<class Array>
void save(Array const& a, std::string const& name, format f) {
	std::ofstream ofs(name);
	*[&]() -> std::unique_ptr<barch::polymorphic_oarchive> {
		switch(f) {
		case xml: return std::make_unique<barch::polymorphic_xml_oarchive>(ofs);
		case txt: return std::make_unique<barch::polymorphic_text_oarchive>(ofs);
		case bin: return std::make_unique<barch::polymorphic_binary_oarchive>(ofs);
		}
		BOOST_MULTI_UNSWITCH();
	}() << bs11n::make_nvp("root", a);
	assert(ofs);
}

template<class Array>
void save(Array const& a, std::experimental::filesystem::path p) {
	if(p.extension() == ".xml")
		return save(a, p.string(), xml);
	else if(p.extension() == ".txt")
		return save(a, p.string(), txt);
	else
		return save(a, p.string(), bin);
}

template<class Array>
void load(Array& a, std::string const& name, format f) {
	std::ifstream ifs(name);
	*[&]() -> std::unique_ptr<barch::polymorphic_iarchive> {switch(f){
		case xml: return std::make_unique<barch::polymorphic_xml_iarchive   >(ifs);
		case txt: return std::make_unique<barch::polymorphic_text_iarchive  >(ifs);
		case bin: return std::make_unique<barch::polymorphic_binary_iarchive>(ifs);
	}BOOST_MULTI_UNSWITCH(); }() >> bs11n::make_nvp("root", a);
	assert(ifs);
}

template<class Array>
void save_xml(Array const& a, std::string const& name) {
	std::ofstream ofs(name);
	barch::xml_oarchive(ofs) << bs11n::make_nvp("root", a);
}

namespace multi = boost::multi;

int main() {
	multi::array<double, 2> const arrD2d = {
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
		{7.0, 8.0, 9.0},
	};
	save(arrD2d, "arrD2d.xml");

	multi::array<double, 2> arrD2d_copy;
	load(arrD2d_copy, "arrD2d.xml", format::xml);
	assert(arrD2d_copy == arrD2d);
}

#undef BOOST_MULTO_SWITCH
