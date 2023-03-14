// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa

#ifndef MPI3_OSTREAM_HPP
#define MPI3_OSTREAM_HPP

#include "../mpi3/communicator.hpp"
#include "../mpi3/process.hpp"

#include <boost/icl/split_interval_map.hpp>  // needs Boost.ICL

#include <iomanip>
#include <iostream>
#include <sstream>

#ifdef __clang__
#define CALLABLE_WHEN(ConsumedORUnconsumed)          [[clang::callable_when(ConsumedORUnconsumed)]]  // NOLINT(cppcoreguidelines-macro-usage)
#define CONSUMABLE(ConsumedORUnconsumed)             [[clang::consumable(ConsumedORUnconsumed)]]  // NOLINT(cppcoreguidelines-macro-usage)
#define RETURN_TYPESTATE(ConsumedORUnconsumed)       [[clang::return_typestate(ConsumedORUnconsumed)]]  // NOLINT(cppcoreguidelines-macro-usage)
#define SET_TYPESTATE(ConsumedORUnconsumedORUnknown) [[clang::set_typestate(ConsumedORUnconsumedORUnknown)]]  // NOLINT(cppcoreguidelines-macro-usage)
#else
#define CALLABLE_WHEN(ConsumedORUnconsumed)  // NOLINT(cppcoreguidelines-macro-usage)
#define CONSUMABLE(ConsumedORUnconsumed)  // NOLINT(cppcoreguidelines-macro-usage)
#define RETURN_TYPESTATE(ConsumedORUnconsumed)  // NOLINT(cppcoreguidelines-macro-usage)
#define SET_TYPESTATE(ConsumedORUnconsumedORUnknown)  // NOLINT(cppcoreguidelines-macro-usage)
#endif

namespace boost {
namespace mpi3 {

struct ostream : public std::ostream {  // NOLINT(fuchsia-multiple-inheritance) bug in clang-tidy 12?
	class streambuf : public std::stringbuf {  // http://stackoverflow.com/a/2212940/225186
		communicator& comm_;
		std::ostream& output;
		std::string   msg_;

		bool doing_table      = false;
		bool doing_formatting = false;

	 public:
		explicit streambuf(communicator& comm, std::ostream& strm = std::cout)
		: comm_{comm}, output{strm} {}

		int sync() override {
			// following code can be improved by a custom reduce operation
			if(comm_.at_root()) {
				boost::icl::interval_map<int, std::string> messages;
				messages.insert(std::make_pair(0, str()));
				for(int i = 1; i != comm_.size(); ++i) {  // NOLINT(altera-unroll-loops) TODO(correaa) use algorithm
					match m = comm_.matched_probe(i);
					msg_.resize(static_cast<std::size_t>(m.count<char>()));
					m.receive(msg_.begin());
					messages.insert(std::make_pair(i, msg_));
				}
				if(
					std::all_of(messages.begin(), messages.end(), [](auto const& e) { return e.second.empty(); }) or
					std::all_of(messages.begin(), messages.end(), [](auto const& e) { return e.second.empty() or e.second.back() == '\n'; })
				) {
					if(not doing_formatting) {
						doing_formatting = true;
						output << '\n';
					}
					if(doing_table) {
						doing_table = false;
						output << '\n';
					}
					collapse_lines(messages);
				} else if(std::all_of(messages.begin(), messages.end(), [](auto const& e) { return e.second.empty() or e.second.back() == '\t'; })) {
					if(not doing_formatting) {
						doing_formatting = true;
						output << '\n';
					}
					if(not doing_table) {
						output  // <<"\e[1;31m"
							<< std::setw(18) << (comm_.name() + "↓");
						for(int i = 0; i != comm_.size(); ++i) {  // NOLINT(altera-unroll-loops) TODO(correaa) use algorithm
							output << std::setw(20) << ("↓" + std::to_string(i) + "↓");
						}
						output  //<<"\e[0m"
							<< '\n';
						doing_table = true;
					}
					headed_row(messages);
				} else {
					doing_formatting = false;
					unformatted_one_or_all(messages);
				}
			} else {
				comm_.send_n(str().begin(), str().size(), 0);
			}
			str("");
			output.flush();
			comm_.barrier();
			return 0;
		}

	 private:
		void collapse_lines(boost::icl::interval_map<int, std::string> const& messages) const {
			for(auto const& m : messages) {  // NOLINT(altera-unroll-loops) use algorithm
				std::string range = comm_.name();
				if(static_cast<int>(size(m.first)) < static_cast<int>(comm_.size())) {
					if(size(m.first) == 1) {
						range += ("[" + std::to_string(lower(m.first)) + "]");
					} else {
						range += ("[" + std::to_string(lower(m.first)) + "-" + std::to_string(upper(m.first)) + "]");
					}
				}
				output << /*"\e[1;32m"<<*/ std::setw(16) << std::left << range;
				output << "→ " << /*"\e[0m"<<*/ m.second;
			}
			if(messages.iterative_size() > 1) {
				output << '\n';
			}
		}
		void headed_row(boost::icl::interval_map<int, std::string> const& messages) const {
			std::size_t last_idx2 = 0;
			std::size_t last_idx  = 0;
			while(  // NOLINT(altera-unroll-loops) use algorithm
				std::all_of(messages.begin(), messages.end(), [last_idx2, &messages](auto const& e) { return last_idx2 != e.second.size() and e.second[last_idx2] == messages.begin()->second[last_idx2]; })
			) {
				if(messages.begin()->second[last_idx2] == ' ') {
					last_idx = last_idx2 + 1;
				}
				++last_idx2;
			}
			output << std::setw(16) << std::left << std::setfill(' ') << messages.begin()->second.substr(0, last_idx);

			for(auto const& m : messages) {
				for(auto i = m.first.lower(); i != m.first.upper() + 1; ++i) {  // NOLINT(altera-id-dependent-backward-branch,altera-unroll-loops) use algorithms
					output
						<< std::setw(16) << std::setfill(' ') << std::left
						<< m.second.substr(last_idx, m.second.size() - last_idx - 1)  // .substr(last_idx, last_idx + 5 /*m.second.size() - 4*/)
						;
				}
			}
			output << '\n'
			       << std::flush;
		}
		void unformatted_one_or_all(boost::icl::interval_map<int, std::string> const& messages) const {
			if(messages.iterative_size() == 1) {
				output << messages.begin()->second << '\n';
			} else {
				for(auto const& m : messages) {
					for(auto i = m.first.lower(); i != m.first.upper() + 1; ++i) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch) TODO(correaa) use algorithm
						output << m.second;
					}
				}
			}
			output << std::flush;
		}
	};

	ostream(ostream const&)            = delete;
	ostream& operator=(ostream const&) = delete;

	ostream(ostream&&)            = delete;
	ostream& operator=(ostream&&) = delete;

 private:
	streambuf buffer;

 public:
	explicit ostream(communicator& comm, std::ostream& os = std::cout)
	: std::ostream(&buffer), buffer(comm, os) {}
	~ostream() override { flush(); }

	// template<class T = void>
	// struct CONSUMABLE(unconsumed)
	// iterator {
	// 	using iterator_category = std::output_iterator_tag;
	// 	using value_type = void;
	// 	using difference_type = std::ptrdiff_t;
	// 	using pointer = void;
	// 	using reference = void;
	// 	using ostream_type = ostream;

	// 	RETURN_TYPESTATE(unconsumed)
	// 	explicit iterator(package_oarchive& oa) : out_{&oa} {}

	// 	iterator& operator=(T const& value) { *out_ << value; return *this; }

	// 	RETURN_TYPESTATE(consumed) CALLABLE_WHEN(unconsumed)
	// 	auto operator*() -> iterator& {return *this;}

	// 	RETURN_TYPESTATE(unconsumed)
	// 	auto operator++() -> iterator& {return *this;}
	// 	RETURN_TYPESTATE(unconsumed)
	// 	auto operator++(int) -> iterator& {return *this;}

	//  private:
	//  	ostream_type* out_;
	// };
};

}  // end namespace mpi3
}  // end namespace boost

#endif
